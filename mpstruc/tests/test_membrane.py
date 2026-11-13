from types import SimpleNamespace

from mpstruc.membrane import _apply_isotonic, _isotonic, explain, features, predict


def _facts(**kw):
    f = {"pdb": "9NEW", "released": "2025-01-10", "uniprots": ["P1"], "tm_uniprots": ["P1"], "keywords": ["KW-0812"],
         "detergents": ["LMN"], "lipids": [], "nanodisc": False, "ipr": ["IPR1"], "methods": "ELECTRON MICROSCOPY",
         "title": "Cryo-EM structure of a sodium channel", "n_protein": 1, "n_unmapped": 0}
    f.update(kw)
    return f


def _ctx(**kw):
    c = SimpleNamespace(protein_dates={"P1": ["2020-01-01"]}, protein_example={"P1": "1ABC"},
                        fam_stats={"IPR1": [9, 1]}, tcdb={"P1": ["1.A.1.1.1"]})
    c.__dict__.update(kw)
    return c


def test_isotonic_is_monotone_and_inside_0_1():
    preds = [(i / 100, int(i > 60)) for i in range(100)] + [(0.95, 0), (0.05, 1)]
    knots = _isotonic(preds, bins=10)
    vals = [y for _, y in knots]
    assert vals == sorted(vals) and 0 < vals[0] and vals[-1] < 1
    assert _apply_isotonic(knots, 0.0) == vals[0] and _apply_isotonic(knots, 1.0) == vals[-1]
    assert _apply_isotonic(knots, 0.2) < _apply_isotonic(knots, 0.9)


def test_features_use_only_earlier_mpstruc_entries():
    x, why = features(_facts(), _ctx(), prior=0.1)
    assert x["same_protein"] == 1 and why["same_protein"] == "1ABC"
    x, _ = features(_facts(released="2019-06-01"), _ctx(), prior=0.1)
    assert x["same_protein"] == 0  # mpstruc only listed it later


def test_family_vote_excludes_the_structure_itself():
    x_own, _ = features(_facts(), _ctx(fam_stats={"IPR1": [1, 0]}), prior=0.1, own_label=1)
    assert x_own["family"] == 0  # its own count removed, nothing left
    x_other, _ = features(_facts(), _ctx(fam_stats={"IPR1": [1, 0]}), prior=0.1)
    assert x_other["family"] > 0


def test_tcdb_and_explanation():
    model = {"bias": -3.0, "weights": {"kw_tm": 3.0, "detergent": 1.0, "tcdb_transporter": 1.0, "same_protein": 1.0,
                                       "kw_KW-0964": -2.0}}
    x, why = features(_facts(), _ctx(), prior=0.1)
    assert x["tcdb_transporter"] == 1
    p = predict(x, model)
    reasons = explain(x, why, _facts(), p, model)
    assert p > 0.5 and reasons[0].startswith("UniProt lists")
    assert any("detergent (LMN)" in r for r in reasons) or any("TCDB" in r for r in reasons)
    # A secreted protein with no membrane keyword reads the other way round.
    f = _facts(tm_uniprots=[], keywords=["KW-0964"], detergents=[])
    x, why = features(f, _ctx(tcdb={}), prior=0.1)
    p = predict(x, model)
    reasons = explain(x, why, f, p, model)
    assert p < 0.5 and any("crossing the membrane" in r for r in reasons)


def test_tm_coverage_detects_missing_membrane_segments():
    from mpstruc.membrane import tm_coverage
    seg = {"P05106": [[719, 741]], "P51681": [[31, 58], [69, 89]]}
    ecto = [{"align": {"P05106": [[1, 1, 788]]}, "unobs": [[508, 788]]}]
    assert tm_coverage(ecto, seg) == (1, 0)          # in the construct, but not resolved
    truncated = [{"align": {"P05106": [[1, 1, 600]]}, "unobs": []}]
    assert tm_coverage(truncated, seg) == (0, 0)     # construct stops before the membrane segment
    ccr5 = [{"align": {"P51681": [[4, 2, 222]]}, "unobs": [[1, 20]]}]
    assert tm_coverage(ccr5, seg) == (2, 2)


def _barrel(**kw):
    """An outer-membrane beta-barrel as UniProt describes it: keywords but no mapped segments."""
    f = dict(uniprots=["Q8A1G1"], tm_uniprots=["Q8A1G1"], keywords=["KW-0472", "KW-0812", "KW-0998", "KW-1134"],
             has_tm_protein=False, tm_unmapped=["Q8A1G1"], tm_in_construct=0, tm_modelled=0, detergents=[],
             title="Bacteroides SusCD open with plug")
    return _facts(**{**f, **kw})


def test_beta_barrel_features():
    ctx = _ctx(tcdb={"Q8A1G1": ["1.B.14.6.1", "2.A.1.1.1"]}, protein_dates={}, protein_example={})
    x, why = features(_barrel(), ctx, prior=0.1)
    assert x["kw_strand"] == 1 and x["kw_outer"] == 1 and x["tcdb_barrel"] == 1
    assert x["tm_unmapped"] == 1 and x["tm_missing"] == 0  # no segments to miss: not an ectodomain
    assert why["tcdb"][0].startswith("1.B.")
    # A helical protein with mapped segments is not "unmapped", even if one partner has no segments.
    x, _ = features(_facts(has_tm_protein=True, tm_unmapped=["P2"], tm_in_construct=7, tm_modelled=7), _ctx(), prior=0.1)
    assert x["tm_unmapped"] == 0 and x["kw_strand"] == 0 and x["tm_modelled"] == 1
    # Coverage unknown: the indicator is left out rather than guessed.
    x, _ = features(_facts(), _ctx(), prior=0.1)
    assert "tm_unmapped" not in x


def test_supporting_reason_never_follows_but():
    model = {"bias": -4.0, "weights": {"kw_tm": 1.0, "kw_strand": 0.8, "t:domain": -0.5, "detergent": 1.0}}
    f = _barrel(title="Plug domain")
    x, why = features(f, _ctx(tcdb={}, protein_dates={}, fam_stats={}), prior=0.1)
    p = predict(x, model)
    reasons = explain(x, why, f, p, model)
    assert p < 0.5
    assert not any(r.startswith("but ") for r in reasons)
    assert reasons[0].startswith("the title mentions 'domain'") and ", although UniProt lists" in reasons[0]
    # Nothing against it at all: still no "but".
    model["weights"].pop("t:domain")
    p = predict(x, model)
    reasons = explain(x, why, f, p, model)
    assert p < 0.5 and reasons[0].startswith("the evidence for it is limited: UniProt")
    # A high score with evidence against it keeps "but" for that evidence only.
    model = {"bias": -1.0, "weights": {"kw_tm": 3.0, "detergent": 1.0, "t:domain": -0.5}}
    p = predict(x, model)
    reasons = explain(x, why, f, p, model)
    assert p > 0.5 and reasons[-1] == "but the title mentions 'domain'"


def test_group_recall_by_kind():
    from mpstruc.membrane import group_recall
    preds = [("1AAA", 0.9, 1, {"keywords": ["KW-0998"]}), ("1BBB", 0.2, 1, {"keywords": ["KW-0998"]}),
             ("1CCC", 0.7, 1, {"keywords": []}), ("1DDD", 0.9, 0, {"keywords": ["KW-0998"]})]
    groups = {"1AAA": "TRANSMEMBRANE PROTEINS: BETA-BARREL", "1BBB": "TRANSMEMBRANE PROTEINS: BETA-BARREL",
              "1CCC": "TRANSMEMBRANE PROTEINS: ALPHA-HELICAL"}
    r = group_recall(preds, groups)
    assert r["TRANSMEMBRANE PROTEINS: BETA-BARREL"] == {"n": 2, "recall_pct": 50.0}
    assert r["TRANSMEMBRANE PROTEINS: ALPHA-HELICAL"] == {"n": 1, "recall_pct": 100.0}
    assert r["outer membrane (UniProt)"] == {"n": 2, "recall_pct": 50.0}  # non-membrane rows are not counted


def test_barrel_domain_counts_only_when_resolved():
    from mpstruc.membrane import barrel_domains, raw_facts
    named = [[116, 238, "TBDR plug"], [243, 1003, "TBDR beta-barrel"], [210, 338, "OmpA-like"]]
    assert barrel_domains(named) == [[243, 1003]]  # the plug and OmpA-like domains stay outside the membrane
    feats = {"Q8A1G1": {"tm": [], "named": named[:2]}}
    det = {"polymer_entities": [{"entity_poly": {"rcsb_entity_polymer_type": "Protein"},
                                 "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["Q8A1G1"]}}]}
    ext = {"kw": {"Q8A1G1": ["KW-0812", "KW-1134", "KW-0998"]}, "released": "2026-11-04"}
    whole = [{"align": {"Q8A1G1": [[1, 105, 899]]}, "unobs": [[493, 495], [540, 577]]}]
    f = raw_facts("38FI", ext, {}, det, whole, feats)
    assert f["has_barrel"] and f["barrel_modelled"] == 1 and f["tm_unmapped"] == []
    x, _ = features(f, _ctx(tcdb={}, protein_dates={}), prior=0.1)
    assert x["tm_modelled"] == 1 and x["barrel_modelled"] == 1 and x["tm_missing"] == 0 and x["tm_unmapped"] == 0
    plug_only = [{"align": {"Q8A1G1": [[1, 105, 140]]}, "unobs": []}]  # residues 105-244: plug, not the barrel
    f = raw_facts("9PLG", ext, {}, det, plug_only, feats)
    assert f["barrel_in_construct"] == 0 and f["barrel_modelled"] == 0
    x, _ = features(f, _ctx(tcdb={}, protein_dates={}), prior=0.1)
    assert x["tm_modelled"] == 0 and x["tm_missing"] == 1  # has a barrel, but not in this structure
    disordered = [{"align": {"Q8A1G1": [[1, 105, 899]]}, "unobs": [[100, 899]]}]
    assert raw_facts("9DIS", ext, {}, det, disordered, feats)["barrel_modelled"] == 0


def _one_protein(acc, kws, title="x"):
    det = {"polymer_entities": [{"entity_poly": {"rcsb_entity_polymer_type": "Protein"},
                                 "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": [acc]}}],
           "struct": {"title": title}}
    return det, {"kw": {acc: kws}, "released": "2025-01-10"}


def test_explain_keeps_reason_limit_with_many_resolved_segments():
    model = {"bias": 0.0, "weights": {"kw_tm": 1.0, "tm_modelled": 2.0, "detergent": 1.0, "same_protein": 0.5,
                                      "family": 0.5, "lipid": 0.4, "nanodisc": 0.4, "t:sodium": -0.5, "kw_strand": 0.3}}
    f = _facts(tm_modelled=16, tm_in_construct=16, has_tm_protein=True, lipids=["PCW"], nanodisc=True,
               keywords=["KW-0812", "KW-1134"], methods="X-RAY DIFFRACTION")
    x, why = features(f, _ctx(), prior=0.1)
    p = predict(x, model)
    reasons = explain(x, why, f, p, model)
    assert p >= 0.5 and len(reasons) == 4 and reasons[-1].startswith("but ")
    assert any("16 membrane-crossing segments are resolved" in r for r in reasons)
    assert len(explain(x, why, f, p, model, n=2)) == 3


def test_aerolysin_prepore_region_is_not_a_membrane_barrel():
    from mpstruc.membrane import barrel_domains, raw_facts
    note = ("Part of the transmembrane beta-barrel after proteolytic activation of the toxin and "
            "insertion into the host membrane")
    named = [[256, 288, note]]
    assert barrel_domains(named, ["KW-0812", "KW-0964"]) == []
    det, ext = _one_protein("P09167", ["KW-0812", "KW-0964", "KW-0800"], "Aerolysin Y221G - prepore")
    cov = [{"align": {"P09167": [[1, 1, 470]]}, "unobs": []}]
    f = raw_facts("9FMX", ext, {}, det, cov, {"P09167": {"tm": [], "named": named}})
    assert f["barrel_modelled"] == 0 and f["tm_unmapped"] == []  # secreted toxin: no "could not be checked" lean
    x, why = features(f, _ctx(tcdb={}, protein_dates={}), prior=0.1)
    model = {"bias": 0.0, "weights": {"tm_modelled": 2.0, "barrel_modelled": 0.5, "tm_unmapped": 0.3}}
    reasons = explain(x, why, f, predict(x, model), model)
    assert not any("beta-barrel that crosses the membrane" in r for r in reasons)
    assert x["tm_unmapped"] == 0


def test_unreviewed_porin_without_transmembrane_keyword_counts_as_barrel():
    from mpstruc.membrane import raw_facts
    det, ext = _one_protein("L0RVU0", ["KW-0002", "KW-0472", "KW-0479", "KW-0732"], "chitoporin")
    cov = [{"align": {"L0RVU0": [[1, 26, 325]]}, "unobs": []}]
    f = raw_facts("8XWO", ext, {}, det, cov, {"L0RVU0": {"tm": [], "named": [[12, 344, "Porin"]]}})
    assert f["has_barrel"] and f["barrel_modelled"] == 1


def test_barrel_after_tonb_plug_when_uniprot_names_only_the_plug():
    from mpstruc.membrane import BARREL_AFTER_PLUG, barrel_domains, raw_facts
    named = [[146, 252, "TonB-dependent receptor plug"]]
    assert barrel_domains(named, ["KW-0998", "KW-1134"]) == [[253, 252 + BARREL_AFTER_PLUG]]
    assert barrel_domains(named, ["KW-0472"]) == []  # not marked outer-membrane: no guess
    assert barrel_domains(named + [[300, 900, "TBDR beta-barrel"]], ["KW-0998"]) == [[300, 900]]
    det, ext = _one_protein("Q64UE4", ["KW-0812", "KW-0998", "KW-1134"])
    feats = {"Q64UE4": {"tm": [], "named": named}}
    f = raw_facts("8YPT", ext, {}, det, [{"align": {"Q64UE4": [[1, 115, 872]]}, "unobs": []}], feats)
    assert f["barrel_modelled"] == 1 and f["tm_unmapped"] == []
    f = raw_facts("9PLG", ext, {}, det, [{"align": {"Q64UE4": [[1, 115, 140]]}, "unobs": []}], feats)
    assert f["barrel_modelled"] == 0 and f["has_barrel"]  # only the plug is in the structure


def test_tm_unmapped_only_for_outer_membrane_proteins():
    from mpstruc.membrane import raw_facts
    feats = {"P1": {"tm": [], "named": []}}
    cov = [{"align": {"P1": [[1, 1, 500]]}, "unobs": []}]
    for kws, expect in [(["KW-0812", "KW-0998"], ["P1"]), (["KW-0812", "KW-1134"], ["P1"]),
                        (["KW-0812"], []),                              # e.g. latrotoxin, NTNH
                        (["KW-0812", "KW-1134", "KW-0964"], [])]:       # secreted pore-forming toxin
        det, ext = _one_protein("P1", kws)
        assert raw_facts("9ABC", ext, {}, det, cov, feats)["tm_unmapped"] == expect, kws


def test_family_evidence_comes_from_the_membrane_protein():
    from mpstruc.membrane import raw_facts
    det = {"polymer_entities": [{"entity_poly": {"rcsb_entity_polymer_type": "Protein"},
                                 "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": [a]}} for a in ("BAR", "UBL")]}
    ext = {"kw": {"BAR": ["KW-0812", "KW-0998", "KW-1134"], "UBL": []}, "released": "2025-01-10"}
    fam = {"ipr": {"BAR": ["IPR_BARREL"], "UBL": ["IPR_UBL"]}}
    f = raw_facts("8YPT", ext, fam, det, None, None)
    assert f["ipr"] == ["IPR_BARREL", "IPR_UBL"] and f["ipr_membrane"] == ["IPR_BARREL"]
    x, why = features(f, _ctx(fam_stats={"IPR_UBL": [0, 396], "IPR_BARREL": [3, 1]}, tcdb={}), prior=0.1)
    assert why["family"][0] == "IPR_BARREL" and x["family"] > 0
    f["ipr_membrane"] = []  # no membrane protein in it: any family may speak
    _, why = features(f, _ctx(fam_stats={"IPR_UBL": [0, 396], "IPR_BARREL": [3, 1]}, tcdb={}), prior=0.1)
    assert why["family"][0] == "IPR_UBL"


def test_one_tcdb_sentence_for_outer_membrane_barrels():
    model = {"bias": -6.0, "weights": {"tcdb_transporter": 0.4, "tcdb_barrel": 0.2, "m:xray": -1.0}}
    f = _barrel(methods="X-RAY DIFFRACTION")
    x, why = features(f, _ctx(tcdb={"Q8A1G1": ["1.B.55.3.1"]}, protein_dates={}, fam_stats={}), prior=0.1)
    p = predict(x, model)
    reasons = explain(x, why, f, p, model, n=5)
    tc = [r for r in reasons if "TCDB" in r]
    assert len(tc) == 1 and "outer-membrane barrel protein (TC 1.B.55.3.1)" in tc[0]


def test_false_positive_split():
    from mpstruc.membrane import false_positive_split
    rows = [("8RM3", 0.65, 0, {"keywords": ["KW-0998", "KW-1134", "KW-0812"]}),
            ("9ARJ", 0.52, 0, {"keywords": ["KW-0812", "KW-0964"]}),
            ("1AAA", 0.9, 1, {"keywords": ["KW-0998", "KW-1134"]}), ("1BBB", 0.2, 0, {"keywords": []})]
    assert false_positive_split(rows) == {"total": 2, "outer_membrane_barrel": 1, "other": 1,
                                          "outer_membrane_barrel_ids": ["8RM3"]}


# ---------------------------------------------------------------- the curators' rule

def _complex(tm_len=21, n_tm=7, tm_unobs=(), partner_len=400, tm_kw=("KW-0812",)):
    """A membrane receptor (P1, n_tm segments of tm_len residues) bound to a big secreted partner (P2)."""
    from mpstruc.membrane import core_tm, tm_by_protein
    segs = [[20 + 40 * i, 20 + 40 * i + tm_len - 1] for i in range(n_tm)]
    cov = [{"align": {"P1": [[1, 1, 300]]}, "unobs": [list(u) for u in tm_unobs]},
           {"align": {"P2": [[1, 1, partner_len]]}, "unobs": []}]
    by, total = tm_by_protein(cov, {"P1": segs, "P2": []})
    core = core_tm(by, total, {"P1"} if tm_kw else set())
    modelled = len(by["P1"]["segments"])
    return _facts(uniprots=["P1", "P2"], tm_uniprots=["P1"] if tm_kw else [], keywords=["KW-0812", "KW-0964"],
                  detergents=[], has_tm_protein=True, tm_in_construct=n_tm, tm_modelled=modelled, barrel_modelled=0,
                  core_tm=core, title="Protein Complex(local refinement)")


# Family and 'secreted' evidence strong enough to pull the model alone well below 50%.
RULE_MODEL = {"bias": -2.0, "weights": {"kw_tm": 1.0, "tm_modelled": 2.0, "family": 1.0,
                                        "kw_KW-0964": -3.0}, "tm_floor": 0.7}


def test_tm_by_protein_counts_resolved_segments_and_share():
    from mpstruc.membrane import core_tm, tm_by_protein
    cov = [{"align": {"P1": [[1, 1, 100]]}, "unobs": [[60, 100]]}, {"align": {"P2": [[1, 1, 300]]}, "unobs": [[1, 100]]}]
    by, total = tm_by_protein(cov, {"P1": [[10, 30], [40, 52], [70, 90]]})
    assert by["P1"] == {"segments": [[10, 30], [40, 52]], "residues": 59}  # the third is not modelled
    assert by["P2"]["residues"] == 200 and total == 259
    core = core_tm(by, total, {"P1"})
    # P1 is 23% of the residues modelled: a minor partner, so it does not carry the rule
    # (and only one of its two resolved segments is long enough).
    assert core["protein"] is None and core["partner_segments"] == 1
    assert core_tm(by, 200, {"P1"})["protein"] == "P1"  # as a main protein it would


def test_floor_applies_to_resolved_segments_of_a_main_protein():
    from mpstruc.membrane import predict_model
    f = _complex()
    assert f["core_tm"]["protein"] == "P1" and f["core_tm"]["segments"] == 7
    x, why = features(f, _ctx(fam_stats={"IPR1": [0, 49]}, tcdb={}), prior=0.1)
    assert x["rule:tm_resolved"] == 1 and "tm_core" not in x  # a floor only, not a learned weight
    assert predict_model(x, RULE_MODEL) < 0.5  # the family and 'secreted' pull the model alone down
    p = predict(x, RULE_MODEL)
    assert p == 0.7
    reasons = explain(x, why, f, p, RULE_MODEL)
    # The resolved segments come first; the counter-evidence is still mentioned after "but".
    assert reasons[0].startswith("7 membrane-crossing segments are resolved in the structure (in UniProt P1")
    assert reasons[-1].startswith("but ")
    # A model without a floor (older membrane_model.json) is unchanged.
    assert predict(x, {k: v for k, v in RULE_MODEL.items() if k != "tm_floor"}) == predict_model(x, RULE_MODEL)


def test_floor_does_not_apply_without_a_significant_resolved_part():
    ctx = _ctx(fam_stats={"IPR1": [0, 49]}, tcdb={})
    for f in (_complex(tm_len=12),                  # short segments: signal peptides, odd annotations
              _complex(n_tm=2),                     # one or two segments: pili, fibrils, anchors
              _complex(tm_unobs=[(1, 300)]),        # not resolved: only the outside part is in the structure
              _complex(partner_len=2000),           # in a minor partner of a big soluble complex
              _complex(tm_kw=())):                  # UniProt does not mark the protein transmembrane
        x, _ = features(f, ctx, prior=0.1)
        assert x["rule:tm_resolved"] == 0
        assert predict(x, RULE_MODEL) < 0.5


def test_resolved_barrel_counts_for_the_rule():
    f = _barrel(has_barrel=True, barrel_modelled=1, tm_modelled=0, tm_unmapped=[],
                core_tm={"protein": None, "segments": 0, "share": 0.0, "strands": 0, "partner_segments": 0})
    x, why = features(f, _ctx(tcdb={}, protein_dates={}), prior=0.1)
    assert x["rule:tm_resolved"] == 1 and why["rule"]["barrel"]
    model = {"bias": -6.0, "weights": {"tm_modelled": 1.0}, "tm_floor": 0.7}
    p = predict(x, model)
    assert p == 0.7 and explain(x, why, f, p, model)[0].startswith("the beta-barrel that crosses the membrane is resolved")


def test_fit_does_not_learn_the_rule_indicator():
    from mpstruc.membrane import fit
    X = [{"rule:tm_resolved": float(i % 2), "a": float(i % 2)} for i in range(40)]
    m = fit(X, [i % 2 for i in range(40)], steps=50)
    assert "rule:tm_resolved" not in m["weights"] and "a" in m["weights"]


def test_segments_are_counted_per_chain_not_summed_over_a_polyprotein():
    """A viral polyprotein (one UniProt entry) whose envelope proteins are separate chains,
    each with two anchors, does not have three segments in one protein."""
    from mpstruc.membrane import core_tm, tm_by_protein
    segs = [[110, 130], [140, 160], [510, 530], [540, 560]]  # prM anchors, then E anchors
    cov = [{"align": {"P1": [[1, 100, 70]]}, "unobs": []},    # prM: residues 100-169
           {"align": {"P1": [[1, 400, 170]]}, "unobs": []}]   # E: residues 400-569
    by, total = tm_by_protein(cov, {"P1": segs})
    assert len(by["P1"]["segments"]) == 2
    core = core_tm(by, total, {"P1"})
    assert core["segments"] == 2
    f = _facts(uniprots=["P1"], tm_uniprots=["P1"], keywords=["KW-0812"], detergents=[], has_tm_protein=True,
               tm_in_construct=4, tm_modelled=4, barrel_modelled=0, core_tm=core, title="virus-like particle")
    x, _ = features(f, _ctx(tcdb={}), prior=0.1)
    assert x["rule:tm_resolved"] == 0
    # The same four segments in one chain do count.
    by, total = tm_by_protein([{"align": {"P1": [[1, 100, 470]]}, "unobs": []}], {"P1": segs})
    assert core_tm(by, total, {"P1"})["segments"] == 4


def test_rule_wording_counts_all_segments_of_the_core_protein():
    from mpstruc.membrane import core_tm, tm_by_protein
    segs = [[20, 40], [50, 70], [80, 100], [110, 121]]  # the last is shorter than CORE_MIN_LENGTH
    cov = [{"align": {"P1": [[1, 1, 300]]}, "unobs": []}]
    by, total = tm_by_protein(cov, {"P1": segs})
    core = core_tm(by, total, {"P1"})
    assert core["segments"] == 3 and core["resolved"] == 4
    f = _facts(uniprots=["P1"], tm_uniprots=["P1"], keywords=["KW-0812"], detergents=[], has_tm_protein=True,
               tm_in_construct=4, tm_modelled=4, barrel_modelled=0, core_tm=core, title="a transporter")
    x, why = features(f, _ctx(tcdb={}), prior=0.1)
    p = predict(x, RULE_MODEL)
    # One UniProt protein: no "N of them in ..." and no "one of its main proteins".
    assert explain(x, why, f, p, RULE_MODEL)[0] == "4 membrane-crossing segments are resolved in the structure"
    # With a partner, all four are in P1, so no "3 of them".
    f2 = dict(f, uniprots=["P1", "P2"])
    x, why = features(f2, _ctx(tcdb={}), prior=0.1)
    assert explain(x, why, f2, predict(x, RULE_MODEL), RULE_MODEL)[0] == \
        "4 membrane-crossing segments are resolved in the structure (in UniProt P1, one of its main proteins)"
