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


# --- per-protein evidence table (tm_table): which protein a curator should check ---

def _entity(desc, ups, kind="Protein"):
    return {"entity_poly": {"rcsb_entity_polymer_type": kind}, "rcsb_polymer_entity": {"pdbx_description": desc},
            "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ups}}


LAT_SEGS = [[37, 57], [69, 89], [107, 127], [133, 153], [160, 180], [186, 206],
            [222, 242], [259, 279], [304, 324], [383, 403], [416, 436], [441, 461]]


def _lat_complex():
    """Like 9TIB: a 12-helix transporter fully modelled, and a partner whose single
    membrane helix (84-104) lies before the modelled part (168-628)."""
    det = {"polymer_entities": [_entity("Y+L amino acid transporter 1", ["Q9UM01"]),
                                _entity("4F2 cell-surface antigen heavy chain", ["P08195"])]}
    cov = [{"align": {"P08195": [[1, 1, 630]]}, "unobs": [[1, 167], [235, 237], [629, 630]]},
           {"align": {"Q9UM01": [[1, 1, 511]]}, "unobs": [[1, 30], [469, 469], [490, 511]]}]
    feats = {"Q9UM01": {"tm": LAT_SEGS, "named": []}, "P08195": {"tm": [[84, 104]], "named": []}}
    return det, cov, feats


def test_tm_table_points_to_the_protein_whose_segments_are_resolved():
    from mpstruc.membrane import tm_by_protein, tm_table
    det, cov, feats = _lat_complex()
    proteins = {"Q9UM01": {"name": "Y+L amino acid transporter 1", "seq": "A" * 511},
                "P08195": {"name": "Amino acid transporter heavy chain SLC3A2", "seq": "A" * 630}}
    rows = tm_table(det, cov, feats, chains={"Q9UM01": 1, "P08195": 1}, proteins=proteins)
    assert [r["uniprot"] for r in rows] == ["Q9UM01", "P08195"]  # the one to look at comes first
    lat, hc = rows
    assert lat["name"] == "Y+L amino acid transporter 1" and lat["chains"] == 1 and lat["length"] == 511
    assert "pdbName" not in lat  # the PDB calls it the same
    # Rows are named as UniProt names the protein; what the PDB calls the chain is kept beside it.
    assert (hc["name"], hc["pdbName"]) == ("Amino acid transporter heavy chain SLC3A2", "4F2 cell-surface antigen heavy chain")
    assert lat["modelled"] == [[31, 468], [470, 489]]
    assert (lat["segments"], lat["resolved"], lat["status"]) == (12, 12, "resolved")
    assert lat["segmentRanges"][0] == [37, 57, True] and len(lat["segmentRanges"]) == 12
    assert lat["note"] == "all 12 membrane-crossing segments are in the structure"
    assert (hc["segments"], hc["resolved"], hc["status"]) == (1, 0, "unresolved")
    assert hc["modelled"] == [[168, 234], [238, 628]] and hc["segmentRanges"] == [[84, 104, False]]
    assert hc["note"] == ("its 1 membrane-crossing segment (84-104) is not in the structure: "
                          "only residues 168-628 were modelled")
    # The same criterion as the score: the table's resolved counts are tm_by_protein's.
    by, _ = tm_by_protein(cov, {a: v["tm"] for a, v in feats.items()})
    assert {r["uniprot"]: r["resolved"] for r in rows} == {a: len(d["segments"]) for a, d in by.items()}


def _protein(n, seed):
    import random
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n))


def test_shifted_pdb_numbering_is_lined_up_by_sequence():
    """9TIB CD98hc: the PDB numbers the chain by a 630-residue form (one region 1->1), UniProt's
    entry has 529 residues, so its helix 84-104 is chain residues 185-205, which are modelled."""
    from mpstruc.membrane import check_alignment, checked_coverage, tm_by_protein, tm_coverage, tm_table
    up = _protein(529, 1)
    chain = _protein(101, 2) + up
    assert check_alignment([[1, 1, 630]], chain, up) == ([[102, 1, 529]], "sequence")
    det, _, feats = _lat_complex()
    item = [{"align": {"P08195": [[1, 1, 630]]}, "seq": chain, "chains": [[[1, 167], [235, 237], [629, 630]]]}]
    proteins = {"P08195": {"name": "Amino acid transporter heavy chain SLC3A2", "seq": up}}
    cov = checked_coverage(item, feats, {}, proteins)
    assert cov[0]["align"] == {"P08195": [[102, 1, 529]]} and cov[0]["realigned"] == ["P08195"]
    hc = next(r for r in tm_table(det, cov, feats, proteins=proteins) if r["uniprot"] == "P08195")
    assert (hc["status"], hc["resolved"], hc["length"], hc["realigned"]) == ("resolved", 1, 529, True)
    assert hc["modelled"] == [[67, 133], [137, 527]] and hc["segmentRanges"] == [[84, 104, True]]
    assert hc["note"].startswith("its 1 membrane-crossing segment (84-104) is in the structure (positions lined up by sequence")
    # The score reads the same coverage: the helix counts there too.
    tm = {a: v["tm"] for a, v in feats.items()}
    assert tm_coverage(cov, tm) == (1, 1) and tm_by_protein(cov, tm)[0]["P08195"]["segments"] == [[84, 104]]
    # Unchecked, the PDB's numbering puts the helix before the modelled part.
    assert tm_coverage([{"align": item[0]["align"], "unobs": item[0]["chains"][0]}], tm) == (1, 0)


def test_check_alignment_deletions_missing_parts_homologs_and_nonsense():
    from mpstruc.membrane import check_alignment, seq_regions
    up = _protein(400, 3)
    # A 30-residue loop deleted from the construct, given by the PDB as one gap-free region.
    cut = up[:150] + up[180:]
    assert check_alignment([[1, 1, 370]], cut, up) == ([[1, 1, 150], [151, 181, 220]], "sequence")
    # Correct regions that leave out the part after the deletion: it is added.
    assert check_alignment([[1, 1, 150]], cut, up) == ([[1, 1, 150], [151, 181, 220]], "sequence")
    # Correct and complete regions are used as given, also with a point mutation in the chain.
    mut = up[:99] + ("A" if up[99] != "A" else "C") + up[100:]
    assert check_alignment([[1, 1, 400]], mut, up) == ([[1, 1, 400]], "pdb")
    assert seq_regions(mut, up) == [[1, 1, 400]]  # substitutions are bridged
    # Numbering past the end of the UniProt sequence is cut off.
    assert check_alignment([[1, 11, 400]], up[10:] + "GSGSGSGSGS", up) == ([[1, 11, 390]], "pdb")
    # A homolog (every 5th residue different) the PDB aligned: kept as given.
    hom = "".join(("A" if c != "A" else "C") if i % 5 == 0 else c for i, c in enumerate(up))
    assert check_alignment([[1, 1, 400]], hom, up) == ([[1, 1, 400]], "pdb-homolog")
    # Two unrelated sequences: nothing may be read off the chain.
    assert check_alignment([[1, 1, 400]], _protein(400, 4), up) == ([], "unchecked")
    # A chimera of two related proteins: the PDB's division between them is kept, although
    # the second would also line up with the first one's part; only unclaimed stretches are added.
    other = "".join(("A" if c != "A" else "C") if i % 20 == 0 else c for i, c in enumerate(up))
    chimera = up[:200] + other[200:]
    assert check_alignment([[201, 201, 200]], chimera, other, claimed=set(range(1, 401))) == ([[201, 201, 200]], "pdb")
    assert check_alignment([[1, 1, 100]], up[:300], up, claimed=set(range(1, 101)) | set(range(250, 301))) == \
        ([[1, 1, 100], [101, 101, 149]], "sequence")


def test_checked_coverage_takes_the_best_modelled_copy_of_a_chain():
    """11UF: two of four copies of the fusion chain have the TARP modelled, the first one has not."""
    from mpstruc.membrane import checked_coverage, tm_table
    up = _protein(300, 5)
    feats = {"P1": {"tm": [[20, 40], [120, 140], [220, 240]], "named": []}}
    item = [{"align": {"P1": [[1, 1, 300]]}, "seq": up, "unobs": [[100, 300]],
             "chains": [[[100, 300]], [[1, 5]], [[100, 300]], [[1, 6]]]}]
    cov = checked_coverage(item, feats, {}, {"P1": {"name": "pump", "seq": up}})
    assert cov[0]["unobs"] == [[1, 5]]
    (row,) = tm_table({"polymer_entities": [_entity("pump", ["P1"])]}, cov, feats)
    assert (row["status"], row["resolved"], row["modelled"]) == ("resolved", 3, [[6, 300]])
    # A protein the entity links but the PDB gives no regions for is lined up from the sequences.
    cov = checked_coverage([{"align": {"P1": []}, "seq": "MGS" + up[10:], "chains": [[]]}], feats, {}, {"P1": {"name": "", "seq": up}})
    assert cov[0]["align"] == {"P1": [[4, 11, 290]]}
    # No UniProt sequence to check against: the PDB's regions stand; with none, nothing is known.
    assert checked_coverage(item, feats)[0]["align"] == {"P1": [[1, 1, 300]]}
    cov = checked_coverage([{"align": {"P1": [[1, 1, 300]]}, "seq": _protein(300, 6), "chains": [[]]}], feats, {},
                           {"P1": {"name": "pump", "seq": up}})
    assert cov[0] == {"align": {}, "unobs": [], "unchecked": ["P1"], "realigned": []}
    (row,) = tm_table({"polymer_entities": [_entity("pump", ["P1"])]}, cov, feats, proteins={"P1": {"name": "pump", "seq": up}})
    assert (row["status"], row["modelled"], row["segmentRanges"][0]) == ("unknown", [], [20, 40, None])
    assert row["note"] == ("UniProt maps 3 membrane-crossing segments (20-40, 120-140, 220-240), but this chain's sequence "
                           "could not be lined up with the UniProt entry, so they could not be checked automatically")
    assert checked_coverage(None) is None


def test_tm_table_deleted_uniprot_entry_is_not_called_soluble():
    from mpstruc.membrane import tm_table
    det = {"polymer_entities": [_entity("ATP-binding cassette, sub-family B (MDR/TAP), member 4", ["E7F1E3"])]}
    cov = [{"align": {"E7F1E3": [[1, 1, 1274]]}, "unobs": []}]
    feats = {"E7F1E3": {"tm": [], "named": []}}  # what uniprot.membrane_features answers for an entry that is gone
    (row,) = tm_table(det, cov, feats, proteins={"E7F1E3": None})
    assert (row["status"], row["gone"], row["length"], row["modelled"]) == ("unknown", True, 0, [])
    assert row["name"].startswith("ATP-binding cassette") and "no longer has this entry" in row["note"]
    (row,) = tm_table(det, cov, feats)  # UniProt not asked: as before
    assert row["status"] == "none" and "gone" not in row


def test_tm_table_uses_the_chain_with_most_resolved_segments_and_counts_partial():
    from mpstruc.membrane import tm_table
    det = {"polymer_entities": [_entity("pump fragment", ["P1"]), _entity("pump", ["P1"])]}
    segs = [[10, 30], [40, 60], [70, 90]]
    cov = [{"align": {"P1": [[1, 1, 35]]}, "unobs": []},                 # first chain: one segment
           {"align": {"P1": [[1, 1, 100]]}, "unobs": [[65, 100]]}]       # second: two of three
    (row,) = tm_table(det, cov, {"P1": {"tm": segs, "named": []}}, chains={"P1": 3})
    assert row["chains"] == 3 and row["modelled"] == [[1, 64]] and row["length"] == 0  # UniProt length not known: none made up
    assert (row["segments"], row["resolved"], row["status"]) == (3, 2, "partial")
    assert row["segmentRanges"] == [[10, 30, True], [40, 60, True], [70, 90, False]]
    assert row["note"] == "2 of its 3 membrane-crossing segments are in the structure (not 70-90)"


def test_tm_table_barrel_soluble_and_unlinked_rows_in_order():
    from mpstruc.membrane import TABLE_MAX_SEGMENTS, tm_table
    named = [[116, 238, "TBDR plug"], [243, 1003, "TBDR beta-barrel"]]
    det = {"polymer_entities": [_entity("nanobody", []), _entity("lipoprotein", ["S1"]), _entity("SusC", ["Q8A1G1"]),
                                _entity("many", ["M1"]), _entity("DNA", [], kind="DNA")]}
    many = [[i * 30, i * 30 + 20] for i in range(1, 41)]
    cov = [{"align": {}, "unobs": []}, {"align": {"S1": [[1, 20, 200]]}, "unobs": []},
           {"align": {"Q8A1G1": [[1, 105, 899]]}, "unobs": [[493, 495]]}, {"align": {"M1": [[1, 1, 1300]]}, "unobs": []}]
    feats = {"Q8A1G1": {"tm": [], "named": named}, "S1": {"tm": [], "named": []}, "M1": {"tm": many, "named": []}}
    rows = tm_table(det, cov, feats, kw={"Q8A1G1": ["KW-1134"]})
    assert [(r["uniprot"], r["status"]) for r in rows] == [("M1", "resolved"), ("Q8A1G1", "resolved"),
                                                           ("S1", "none"), (None, "unknown")]
    assert rows[0]["segments"] == 40 and len(rows[0]["segmentRanges"]) == TABLE_MAX_SEGMENTS
    barrel = rows[1]
    assert barrel["barrel"] == {"barrel": "TBDR beta-barrel", "range": [243, 1003], "resolved": True}
    assert barrel["segments"] == 0 and barrel["note"] == "beta-barrel domain 243-1003 is in the structure"
    assert rows[2]["note"] == "no membrane-crossing segments in UniProt" and rows[2]["modelled"] == [[20, 219]]
    assert rows[3]["name"] == "nanobody" and "not linked to UniProt" in rows[3]["note"]  # the DNA is not a row
    # Only the plug of the receptor in the structure: its barrel is not there.
    plug = tm_table(det, [{"align": {"Q8A1G1": [[1, 105, 140]]}, "unobs": []}], feats)
    q = next(r for r in plug if r["uniprot"] == "Q8A1G1")
    assert q["status"] == "unresolved" and q["barrel"]["resolved"] is False
    assert q["note"] == "beta-barrel domain 243-1003 is not in the structure: only residues 105-244 were modelled"
    assert tm_table(None, None, None) == []
    # A barrel that passes the at-least-half rule but is not all there is not called "in the structure" outright.
    part = tm_table(det, [{"align": {"Q8A1G1": [[1, 105, 899]]}, "unobs": [[567, 899]]}], feats)
    q = next(r for r in part if r["uniprot"] == "Q8A1G1")
    assert q["status"] == "resolved" and q["note"] == "part of beta-barrel domain 243-1003 is in the structure (modelled 243-670)"
    # '1 of its 4 ... is', and a protein in a chain that holds two UniProt proteins is marked.
    det2 = {"polymer_entities": [_entity("A,B chimera", ["P1", "P2"])]}
    cov2 = [{"align": {"P1": [[1, 1, 100]], "P2": [[101, 1, 50]]}, "unobs": [[40, 100]]}]
    one, two = tm_table(det2, cov2, {"P1": {"tm": [[10, 30], [45, 65], [70, 90], [95, 99]], "named": []}, "P2": {"tm": [], "named": []}},
                        proteins={"P1": {"name": "Protein A", "seq": "A" * 100}, "P2": {"name": "Protein B", "seq": "A" * 50}})
    assert one["note"] == "1 of its 4 membrane-crossing segments is in the structure (not 45-65, 70-90, 95-99)"
    assert (one["name"], one["fusion"], two["name"], two["fusion"]) == ("Protein A", True, "Protein B", True)
    assert "pdbName" not in one


def test_reason_names_the_proteins_when_there_are_several():
    from mpstruc.membrane import chain_counts, name_proteins, tm_table
    det, cov, feats = _lat_complex()
    table = tm_table(det, cov, feats)
    f = _facts(uniprots=["P08195", "Q9UM01"], tm_uniprots=["P08195", "Q9UM01"], keywords=["KW-0812"], detergents=[],
               has_tm_protein=True, tm_in_construct=13, tm_modelled=12, barrel_modelled=0, title="a complex",
               core_tm={"protein": "Q9UM01", "segments": 12, "resolved": 12, "share": 0.5, "strands": 0, "partner_segments": 0})
    x, why = features(f, _ctx(tcdb={}), prior=0.1)
    p = predict(x, RULE_MODEL)
    assert explain(x, why, f, p, RULE_MODEL, table=table)[0] == \
        "12 membrane-crossing segments are resolved in the structure, all in Y+L amino acid transporter 1 (Q9UM01)"
    # A single protein is not named (nothing to mix up); without a table the wording is as before.
    assert explain(x, why, f, p, RULE_MODEL, table=table[:1])[0].endswith("one of its main proteins)")
    rows = [{"uniprot": a, "name": n, "resolved": k} for a, n, k in
            [("Q8W4S4", "V-type proton ATPase subunit a3", 6), ("P0DH92", "V-type proton ATPase subunit c1", 4),
             ("Q9SZY7", "V-type proton ATPase subunit c''1", 4), ("Q9SZ13", "V-type proton ATPase subunit e2", 2),
             ("P11574", "V-type proton ATPase subunit B1", 0)]]
    assert name_proteins(rows) == "V-type proton ATPase subunit a3 (Q8W4S4) 6, c1 (P0DH92) 4, c''1 (Q9SZY7) 4 and 1 more"
    # Names with a comma (fusion chains as the PDB writes them) are not cut at a shared first word.
    fus = [{"uniprot": "P19491", "name": "Isoform Flip of Glutamate receptor 2", "resolved": 5},
           {"uniprot": "D3ZK93", "name": "Isoform Flip of Glutamate receptor 2,Germ cell-specific gene 1-like protein chimera", "resolved": 4}]
    assert name_proteins(fus) == ("Isoform Flip of Glutamate receptor 2 (P19491) 5, "
                                  "Isoform Flip of Glutamate receptor 2,Germ cel... (D3ZK93) 4")
    # The number said is the sum of the counts listed, not the total over chain entities (14 here).
    assert explain(x, why, dict(f, tm_modelled=14), p, RULE_MODEL, table=fus)[0].startswith(
        "9 membrane-crossing segments are resolved in the structure: Isoform Flip of Glutamate receptor 2 (P19491) 5, ")
    yidc = [{"uniprot": "P25714", "name": "Membrane protein insertase YidC", "resolved": 5}, {"uniprot": "P0AG27", "name": "YibN", "resolved": 0}]
    assert explain(x, why, dict(f, tm_modelled=6), p, RULE_MODEL, table=yidc)[0] == \
        "5 membrane-crossing segments are resolved in the structure, all in Membrane protein insertase YidC (P25714)"
    f18 = dict(f, tm_modelled=16)
    assert explain(x, why, f18, p, RULE_MODEL, table=rows)[0].startswith(
        "16 membrane-crossing segments are resolved in the structure: V-type proton ATPase subunit a3 (Q8W4S4) 6, c1")
    ent = {"polymer_entities": [
        {"rcsb_polymer_entity": {"pdbx_description": "c1"},
         "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["P0DH92"], "auth_asym_ids": list("abcdefghi")}},
        {"rcsb_polymer_entity": {"pdbx_description": "Fab heavy chain"},
         "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": None, "auth_asym_ids": ["H", "I"]}}]}
    assert chain_counts(ent) == {"P0DH92": 9, "Fab heavy chain": 2}


def test_backfill_adds_the_table_and_keeps_everything_else(tmp_path, monkeypatch):
    import json
    from mpstruc import membrane
    path = tmp_path / "candidates.json"
    data = {"generated": "2026-11-13T01:52+00:00", "stats": {"scored": 2}, "candidates": [
        {"pdb": "9TIB", "title": "café", "membrane": 0.945, "membraneWhy": ["old"], "res": 3.5},
        {"pdb": "9XXX", "membrane": 0.3, "membraneWhy": ["kept"]},
        {"pdb": "9GONE", "membrane": 0.2, "membraneWhy": ["x"], "tmEvidence": [{"stale": 1}]}], "inMpstruc": []}
    path.write_text(json.dumps(data, separators=(",", ":")))
    row = [{"uniprot": "Q9UM01", "resolved": 12}]
    monkeypatch.setattr(membrane, "load_model", lambda: {"prior": 0.1})
    monkeypatch.setattr(membrane, "_scored", lambda ids, db, model: iter(
        [("9TIB", 0.9449, ["new"], row), ("9XXX", 0.8, ["moved"], row)]))
    import datetime
    n = membrane.backfill_table(path, today=datetime.date(2026, 11, 13))
    assert n == {"candidates": 3, "with_table": 2, "without_table": 1, "why_changed": 1, "why_kept_score_moved": 1}
    text = path.read_text()
    out = json.loads(text)
    assert text == json.dumps(out, separators=(",", ":"))  # as automation.py writes it
    a, b, c = out["candidates"]
    assert a == {"pdb": "9TIB", "title": "café", "membrane": 0.945, "membraneWhy": ["new"], "res": 3.5, "tmEvidence": row}
    assert b["membraneWhy"] == ["kept"] and b["tmEvidence"] == row  # the score moved: stored score and text stay together
    assert "tmEvidence" not in c
    assert out.pop("tmChecked") == "2026-11-13"  # the day the tables were checked against the PDB and UniProt
    assert {k: v for k, v in out.items() if k != "candidates"} == {k: v for k, v in data.items() if k != "candidates"}
