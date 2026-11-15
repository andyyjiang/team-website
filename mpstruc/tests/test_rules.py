"""The curators' subfamily rules: loading, barrel-strand counting, placement, membrane guard."""
import json
import math
from collections import Counter

import pytest

from mpstruc import barrel, rules
from mpstruc.automation import KB, Known, apply_rules


# ---------------------------------------------------------------- rules file

def test_the_rules_file_has_the_two_curator_rules_verbatim():
    got = {r["subgroup"]: r for r in rules.load_rules()}
    at, omp = got["Outer Membrane Autotransporters"], got["Omp85-TpsB Outer Membrane Transporter Superfamily"]
    assert (at["note"], at["barrelStrands"]) == ("Typically 12-stranded beta-barrel.", 12)
    assert (omp["note"], omp["barrelStrands"]) == ("Typically 16-stranded beta-barrel.", 16)
    assert all(r["addedBy"] and r["addedOn"] for r in got.values())


def test_rules_are_matched_by_subfamily_name_not_key(tmp_path):
    path = tmp_path / "curator_rules.json"
    path.write_text(json.dumps({"rules": [
        {"subgroup": "Outer  Membrane autotransporters", "subgroupKey": "G02-S006", "note": "Typically 12-stranded.", "barrelStrands": 12},
        {"subgroup": "Renamed since", "subgroupKey": "G02-S099", "note": "kept by key", "barrelStrands": 16, "tolerance": 1},
        {"subgroup": "Gone", "subgroupKey": "G09-S001", "note": "no such subfamily"},
        {"subgroup": "Porins", "note": "a rule without a strand count"}]}))
    loaded = rules.load_rules(path)
    names = {"G02-S007": "Outer Membrane Autotransporters", "G02-S099": "Something else", "G02-S003": "Porins"}  # keys shifted
    got = rules.resolve(loaded, names)
    assert set(got) == {"G02-S007", "G02-S099", "G02-S003"}
    assert got["G02-S007"]["barrelStrands"] == 12 and got["G02-S007"]["tolerance"] == 0
    assert rules.unmatched(loaded, names) == ["Gone"]
    assert set(rules.strand_rules(got)) == {"G02-S007", "G02-S099"}
    assert rules.fits(got["G02-S099"], 15) and rules.fits(got["G02-S099"], 17) and not rules.fits(got["G02-S099"], 14)
    assert not rules.fits(got["G02-S007"], 13) and not rules.fits(got["G02-S007"], None)
    site = rules.site_data(got, [{"pdb": "4C00"}], "2026-11-13")
    assert [r["key"] for r in site["rules"]] == ["G02-S003", "G02-S007", "G02-S099"] and site["conflicts"] == [{"pdb": "4C00"}]
    assert rules.load_rules(tmp_path / "missing.json") == []


@pytest.mark.parametrize("bad", [{"rules": [{"subgroup": "X"}]}, {"rules": [{"subgroup": "X", "note": "n", "barrelStrands": "12"}]},
                                 {"rules": [{"subgroup": "X", "note": "n", "barrelStrands": 12, "tolerance": -1}]}, {"rules": {}}])
def test_a_malformed_rule_is_refused(tmp_path, bad):
    path = tmp_path / "curator_rules.json"
    path.write_text(json.dumps(bad))
    with pytest.raises(rules.RuleError):
        rules.load_rules(path)


# ---------------------------------------------------------------- strand counting

def _ring(n, slots=None, m=10, z0=0.0, first=0, phase=0.0):
    """C-alphas of n antiparallel strands of m residues standing side by side (4.8 A apart)
    on a ring with room for `slots` strands; seq numbers leave a gap of 4 between strands."""
    slots = slots or n
    radius = 2.4 / math.sin(math.pi / slots)
    seq, xyz = [], []
    for k in range(n):
        a = phase + 2 * math.pi * k / slots
        for j in range(m):
            z = z0 + 3.3 * (j if (first + k) % 2 == 0 else m - 1 - j)
            seq.append((first + k) * (m + 4) + j + 1)
            xyz.append([radius * math.cos(a), radius * math.sin(a), z])
    return {"entity": "1", "seq": seq, "xyz": xyz}


@pytest.mark.parametrize("n", [8, 12, 16, 22])
def test_counts_the_strands_of_a_closed_barrel(n):
    res = barrel.count_barrels({"chains": {"A": _ring(n)}, "operators": []})
    assert res["chains"]["A"]["count"] == n and res["chains"]["A"]["closed"] is True
    best = barrel.best_barrel(res)
    assert (best["chain"], best["count"], best["method"]) == ("A", n, barrel.METHOD)
    assert best["first"] in (1, 2) and n * 14 - 6 <= best["last"] <= n * 14


def test_an_open_seam_is_still_counted_and_a_flat_sheet_is_not():
    res = barrel.count_barrels({"chains": {"A": _ring(16, slots=17)}, "operators": []})
    assert res["chains"]["A"]["count"] == 16 and res["chains"]["A"]["closed"] is False
    flat = {"entity": "1", "seq": [k * 14 + j + 1 for k in range(10) for j in range(10)],
            "xyz": [[4.8 * k, 0.0, 3.3 * (j if k % 2 == 0 else 9 - j)] for k in range(10) for j in range(10)]}
    assert barrel.count_barrels({"chains": {"A": flat}, "operators": []})["chains"] == {}
    assert barrel.best_barrel({"chains": {}, "oligomers": []}) is None and barrel.best_barrel(None) is None


def test_a_sheet_outside_the_barrel_and_a_helix_are_not_counted():
    ca = _ring(12)
    # A three-stranded sheet 60 A away (a POTRA domain), earlier in the same chain.
    for k in range(3):
        for j in range(8):
            ca["seq"].append(1000 + k * 12 + j)
            ca["xyz"].append([60.0 + 4.8 * k, 0.0, 3.3 * (j if k % 2 == 0 else 7 - j)])
    # ... and an alpha helix (1.5 A rise, 100 degrees per residue).
    for j in range(20):
        ca["seq"].append(2000 + j)
        ca["xyz"].append([-60 + 2.3 * math.cos(math.radians(100 * j)), 2.3 * math.sin(math.radians(100 * j)), 1.5 * j])
    res = barrel.count_barrels({"chains": {"A": ca}, "operators": []})
    assert res["chains"]["A"]["count"] == 12


def test_a_bulge_does_not_split_a_strand():
    ca = _ring(12, m=12)
    # Residue 6 of strand 3 bulges out of the sheet: its run of bridges is broken there.
    i = 3 * 12 + 6
    x, y, z = ca["xyz"][i]
    ca["xyz"][i] = [x * 1.9, y * 1.9, z]
    res = barrel.count_barrels({"chains": {"A": ca}, "operators": []})
    assert res["chains"]["A"]["count"] == 12


def test_a_barrel_made_by_three_chains_is_reported_with_its_chains():
    chains = {c: _ring(4, slots=12, first=0, phase=2 * math.pi * 4 * k / 12) for k, c in enumerate("ABC")}
    res = barrel.count_barrels({"chains": chains, "operators": []})
    assert res["chains"] == {}
    (o,) = res["oligomers"]
    assert (o["count"], o["chains"], o["perChain"]) == (12, ["A", "B", "C"], 4.0)
    best = barrel.best_barrel(res)
    assert best["count"] == 12 and best["chains"] == 3
    # One chain in the deposited model, the trimer from the assembly's operators.
    rot = lambda a: [[[math.cos(a), -math.sin(a), 0.0], [math.sin(a), math.cos(a), 0.0], [0.0, 0.0, 1.0]], [0.0, 0.0, 0.0]]
    one = barrel.count_barrels({"chains": {"A": chains["A"]}, "operators": [rot(0), rot(2 * math.pi / 3), rot(4 * math.pi / 3)]})
    assert one["oligomers"][0]["count"] == 12 and len(one["oligomers"][0]["chains"]) == 3


def test_parse_cif_reads_loops_single_values_and_quotes():
    text = """data_X
#
_entry.id   1ABC
#
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_comp_id
ATOM CA GLY
ATOM 'C A' "AL A"
#
_pdbx_struct_oper_list.id 1
_pdbx_struct_oper_list.name
;multi
line
;
#
"""
    got = barrel.parse_cif(text, ("atom_site", "pdbx_struct_oper_list"))
    assert got["atom_site"] == {"group_PDB": ["ATOM", "ATOM"], "label_atom_id": ["CA", "C A"], "label_comp_id": ["GLY", "AL A"]}
    assert got["pdbx_struct_oper_list"] == {"id": ["1"], "name": ["multi\nline"]}
    assert "entry" not in got


def test_counts_are_cached_and_a_failed_download_is_not(tmp_path, monkeypatch):
    monkeypatch.setattr(barrel.rcsb, "CACHE_DIR", tmp_path)
    calls = []

    def fake(pdb, entities=None):
        calls.append((pdb, entities))
        if pdb == "9ERR":
            raise OSError("network down")
        return {"chains": {"A": _ring(12)}, "operators": []} if pdb == "9AAA" else None

    monkeypatch.setattr(barrel, "fetch_ca", fake)
    got = barrel.barrel_counts({"9AAA": ["1"], "9NON": None, "9ERR": None})
    assert got["9AAA"]["chains"]["A"]["count"] == 12 and got["9NON"] is None and got["9ERR"] is None
    barrel.barrel_counts({"9AAA": ["1"], "9NON": None, "9ERR": None})
    assert [c[0] for c in calls] == ["9AAA", "9NON", "9ERR", "9ERR"]  # only the failed one is asked again


# ---------------------------------------------------------------- placement

def _s(core, uipr=(), citation=()):
    core = set(core)
    return {"uniprot": core, "core": core, "all_uniprot": core, "taxa": {562}, "pfam": set(), "interpro": set(),
            "opm": set(), "uipr": set(uipr), "citation": set(citation), "released": "2015-01-01"}


AT, OMP = "G02-S006", "G02-S007"
RULED = {AT: {"key": AT, "name": "Outer Membrane Autotransporters", "barrelStrands": 12, "tolerance": 1, "note": "Typically 12-stranded beta-barrel."},
         OMP: {"key": OMP, "name": "Omp85-TpsB Outer Membrane Transporter Superfamily", "barrelStrands": 16, "tolerance": 1,
               "note": "Typically 16-stranded beta-barrel."}}


def _kb():
    kb = KB(entries=[])
    rows = [
        Known("4C00", f"{AT}-P0011", "parent", AT, f"{AT}-P0011", f"{AT}-P0011", "2013", _s(["TAMA_ECOLI"], ["IPR000184", "IPR035243"])),
        Known("2QOM", f"{AT}-P0003", "parent", AT, f"{AT}-P0003", f"{AT}-P0003", "2007", _s(["ESPP"], ["IPR005546"])),
        Known("4K3B", f"{OMP}-P0004", "parent", OMP, f"{OMP}-P0004", f"{OMP}-P0004", "2013", _s(["BAMA"], ["IPR000184", "IPR023707"])),
        Known("8US1", f"{OMP}-P0018", "parent", OMP, f"{OMP}-P0018", f"{OMP}-P0018", "2024", _s(["TAMA_PSEAE"], ["IPR000184", "IPR035243"])),
        Known("8US2", f"{OMP}-P0018-R001", "related", OMP, f"{OMP}-P0018", f"{OMP}-P0018", "2024", _s(["TAMA_PSEAE"], ["IPR000184", "IPR035243"])),
    ]
    for k in rows:
        kb.entries.append(k)
        kb.by_pdb.setdefault(k.pdb, []).append(k)
        kb.by_key[k.key] = k
        for u in k.s["core"]:
            kb.by_uniprot.setdefault(u, []).append(k)
        kb.sg_count[k.sg] += 1
    kb.sg_names = {AT: RULED[AT]["name"], OMP: RULED[OMP]["name"], "G03-S001": "Channels"}
    kb.sg_group = {AT: "Beta-barrel transmembrane", OMP: "Beta-barrel transmembrane", "G03-S001": "Alpha-helical transmembrane"}
    return kb


def _cand(pdb, placement, sg, owner=None, rule="partial-protein", uniprot=("TAMA_ECOLI", "TAMB"), membrane=0.95, alts=(AT, OMP)):
    names = _kb().sg_names
    return {"pdb": pdb, "placement": placement, "subgroup": sg, "subgroupName": names.get(sg), "group": "Beta-barrel transmembrane",
            "ownerPdb": owner, "ownerKey": f"{AT}-P0011" if owner == "4C00" else None, "confidence": 0.0, "rule": rule,
            "evidence": ["first reason"], "uniprot": list(uniprot), "membrane": membrane,
            "alternatives": [{"subgroup": k, "name": names[k], "score": 0.5, "evidence": []} for k in alts]}


def _barrel(count, chain="A", **more):
    return {"chain": chain, "count": count, "method": barrel.METHOD, "closed": True, "entity": "1", **more}


def test_sixteen_strands_move_a_tama_complex_from_autotransporters_to_omp85():
    lead = _cand("9TGQ", "review", AT, "4C00")
    follower = _cand("9TGR", "related", AT, "9TGQ", rule="same-paper-batch")
    props = [lead, follower]
    barrels = {"9TGQ": {"TAMA_ECOLI": _barrel(16)}, "9TGR": {"TAMA_ECOLI": _barrel(16, "B")}}
    fams = {"9TGQ": {"ipr": {"TAMA_ECOLI": ["IPR000184", "IPR035243"]}}}
    changes = apply_rules(props, _kb(), RULED, barrels, {"9TGQ": {"TAMA_ECOLI"}}, fams, {})
    assert (lead["subgroup"], lead["subgroupName"], lead["placement"]) == (OMP, RULED[OMP]["name"], "review")
    assert (lead["ownerPdb"], lead["ownerKey"]) == ("8US1", f"{OMP}-P0018")  # same family there, not BamA
    assert lead["rule"] == "barrel-rule" and lead["ruleCheck"]["status"] == "moved" and lead["ruleCheck"]["from"] == AT
    assert lead["barrelStrands"] == {"chain": "A", "count": 16, "method": barrel.METHOD, "closed": True, "uniprot": "TAMA_ECOLI"}
    assert lead["evidence"][0] == "first reason"
    assert ("16 strands counted in the barrel (chain A): Omp85-TpsB Outer Membrane Transporter Superfamily are typically "
            "16-stranded, Outer Membrane Autotransporters typically 12-stranded (curators' rule)") in lead["evidence"][1]
    assert "4C00" in lead["evidence"][2] and "8US1" in lead["evidence"][3]
    assert follower["subgroup"] == OMP and follower["ownerPdb"] == "9TGQ" and follower["barrelStrands"]["chain"] == "B"
    assert {c["pdb"] for c in changes} == {"9TGQ", "9TGR"}
    # Running it again on its own output changes nothing more (the weekly run and the backfill agree).
    again = json.loads(json.dumps(props))
    assert apply_rules(again, _kb(), RULED, barrels, {"9TGQ": {"TAMA_ECOLI"}}, fams, {}) == []
    assert again[0]["subgroup"] == OMP and again[0]["ruleCheck"]["status"] == "fits"


def test_a_count_that_fits_is_noted_and_other_subfamilies_are_left_alone():
    ok = _cand("9BAM", "parent", OMP, rule="new-protein", uniprot=("BAMA2",))
    twelve = _cand("9ATX", "daughter", AT, "2QOM", rule="same-protein", uniprot=("ESPP",))
    fourteen = _cand("9ODD", "parent", AT, rule="new-protein", uniprot=("X",))
    helical = _cand("9HEL", "parent", "G03-S001", rule="new-protein", uniprot=("Y",), alts=("G03-S001",))
    props = [ok, twelve, fourteen, helical]
    barrels = {"9BAM": {"BAMA2": _barrel(15)}, "9ATX": {"ESPP": _barrel(12)}, "9ODD": {"X": _barrel(14)}, "9HEL": {"Y": _barrel(16)}}
    assert apply_rules(props, _kb(), RULED, barrels, {}, {}, {}) == []
    assert ok["ruleCheck"]["status"] == "fits" and ok["placement"] == "parent" and ok["barrelStrands"]["count"] == 15
    assert twelve["ruleCheck"]["status"] == "fits" and twelve["ownerPdb"] == "2QOM"
    assert fourteen["ruleCheck"]["status"] == "differs" and fourteen["subgroup"] == AT  # fits neither rule: only noted
    assert "ruleCheck" not in helical and "barrelStrands" not in helical and helical["evidence"] == ["first reason"]


def test_a_new_parent_stays_a_new_parent_and_an_official_paper_keeps_its_entry():
    parent = _cand("9NEW", "parent", AT, rule="new-protein", uniprot=("TAMA_NEW",))
    related = _cand("9REL", "related", AT, "4C00", rule="same-paper", uniprot=("TAMA_ECOLI",))
    barrels = {"9NEW": {"TAMA_NEW": _barrel(16)}, "9REL": {"TAMA_ECOLI": _barrel(16)}}
    fams = {"9NEW": {"ipr": {"TAMA_NEW": ["IPR000184", "IPR035243"]}}}
    apply_rules([parent, related], _kb(), RULED, barrels, {}, fams, {})
    assert (parent["placement"], parent["subgroup"], parent["ownerPdb"]) == ("parent", OMP, None)
    assert parent["evidence"][-1] == "Closest entry there: 8US1 (same protein family)"
    assert (related["placement"], related["subgroup"], related["ownerPdb"]) == ("related", AT, "4C00")
    assert related["ruleCheck"]["status"] == "differs" and "stays with 4C00" in related["ruleCheck"]["text"]


def test_a_barrel_of_a_substrate_does_not_move_a_bam_complex():
    c = _cand("9HYB", "parent", OMP, rule="new-protein", uniprot=("BAMA2", "ESPP2"))
    apply_rules([c], _kb(), RULED, {"9HYB": {"ESPP2": _barrel(12, "F"), "BAMA2": _barrel(16)}}, {}, {}, {})
    assert c["subgroup"] == OMP and c["ruleCheck"]["status"] == "fits" and c["barrelStrands"]["uniprot"] == "BAMA2"


def test_no_barrel_and_a_low_score_go_to_review_instead_of_a_barrel_subfamily():
    pore = _cand("9AAV", "parent", AT, rule="new-protein", uniprot=("NUP",), membrane=0.171)
    twin = _cand("9AB1", "related", AT, "9AAV", rule="same-paper-batch", uniprot=("NUP",), membrane=0.171)
    potra = _cand("9TRP", "parent", OMP, rule="new-protein", uniprot=("BAMA_BORBU",), membrane=0.305)
    sure = _cand("9SUR", "parent", AT, rule="new-protein", uniprot=("AT2",), membrane=0.9)  # e.g. a named barrel domain, resolved
    small = _cand("9SML", "parent", AT, rule="new-protein", uniprot=("Z",), membrane=0.2)
    props = [pore, twin, potra, sure, small]
    # Every structure was counted ({}: no barrel found); only 9SML has a ring, of five strands.
    barrels = {"9AAV": {}, "9AB1": {}, "9TRP": {}, "9SUR": {}, "9SML": {"Z": _barrel(5)}}
    changes = apply_rules(props, _kb(), RULED, barrels, {"9TRP": {"BAMA_BORBU"}, "9SML": {"Z"}}, {}, {})
    assert (pore["placement"], pore["subgroup"], pore["subgroupName"], pore["confidence"]) == ("review", None, None, 0.0)
    assert pore["rule"] == "barrel-rule-no-barrel" and pore["ruleCheck"]["status"] == "no-barrel"
    assert "none of its proteins is an outer-membrane barrel protein" in pore["evidence"][-1]
    assert twin["subgroup"] is None and twin["placement"] == "related" and twin["ruleCheck"]["status"] == "no-barrel"
    # An outer-membrane protein's soluble part: to review, but the subfamily stays as the suggestion.
    assert (potra["placement"], potra["subgroup"]) == ("review", OMP) and "curator decides" in potra["evidence"][-1]
    assert (sure["placement"], sure["subgroup"]) == ("parent", AT) and "ruleCheck" not in sure
    assert small["placement"] == "review" and "barrelStrands" not in small  # five strands are not a membrane barrel
    assert {c["pdb"] for c in changes} == {"9AAV", "9AB1", "9TRP", "9SML"}


# ---------------------------------------------------------------- membrane score

def _structure(acc="BAMA", kws=("KW-0998",), ipr=("IPR000184",), named=(), tm=()):
    det = {"polymer_entities": [{"entity_poly": {"rcsb_entity_polymer_type": "Protein"},
                                 "rcsb_polymer_entity": {"pdbx_description": "barrel protein"},
                                 "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": [acc]}}],
           "exptl": [{"method": "ELECTRON MICROSCOPY"}], "struct": {"title": "A complex"}}
    ext = {"kw": {acc: list(kws)}, "released": "2026-11-01"}
    fam = {"ipr": {acc: list(ipr)}}
    cov = [{"align": {acc: [[1, 21, 790]]}, "unobs": []}]
    feats = {acc: {"tm": [list(t) for t in tm], "named": [list(n) for n in named]}}
    return ext, fam, det, cov, feats


def test_a_counted_barrel_in_an_outer_membrane_protein_is_a_resolved_barrel():
    from mpstruc.membrane import CORE_MIN_STRANDS, TM_FLOOR, features, predict, raw_facts, tm_table
    ext, fam, det, cov, feats = _structure(named=[(24, 91, "POTRA")])
    counted = {"BAMA": {**_barrel(16), "first": 404, "last": 786}}
    plain = raw_facts("9PXG", ext, fam, det, cov, feats)
    assert plain["barrel_counted"] == {} and not plain["barrel_modelled"] and not plain["has_barrel"]
    f = raw_facts("9PXG", ext, fam, det, cov, feats, counted)
    assert f["barrel_counted"] == {"BAMA": 16} and f["barrel_modelled"] == 1 and f["has_barrel"] is True

    class Ctx:
        tcdb, protein_dates, protein_example, fam_stats = {}, {}, {}, {}
    model = {"bias": -4.0, "weights": {}, "tm_floor": TM_FLOOR}
    x0, _ = features(plain, Ctx, 0.1)
    x1, _ = features(f, Ctx, 0.1)
    assert not x0.get("rule:tm_resolved") and predict(x0, model) < 0.1
    assert x1["rule:tm_resolved"] == 1.0 and x1["barrel_modelled"] == 1.0 and predict(x1, model) == TM_FLOOR
    (row,) = tm_table(det, cov, feats, ext["kw"], counted=counted)
    assert row["status"] == "resolved" and row["barrel"]["strands"] == 16
    assert row["barrel"]["range"] == [424, 806]  # the chain's residues 404-786 in UniProt numbering
    assert row["note"] == "16-strand beta-barrel counted in the structure (UniProt names no barrel domain)"
    # A barrel whose seam is open is a lower bound, and is said to be one.
    (row,) = tm_table(det, cov, feats, ext["kw"], counted={"BAMA": {**_barrel(13, closed=False), "first": 404, "last": 786}})
    assert row["barrel"]["open"] is True and row["barrel"]["strands"] == 13
    assert row["note"] == "beta-barrel of at least 13 strands counted in the structure, its seam open (UniProt names no barrel domain)"
    # Too few strands to cross a membrane as a barrel.
    few = {"BAMA": _barrel(CORE_MIN_STRANDS - 1)}
    assert raw_facts("9PXG", ext, fam, det, cov, feats, few)["barrel_counted"] == {}
    assert tm_table(det, cov, feats, ext["kw"], counted=few)[0]["status"] == "none"


@pytest.mark.parametrize("kws, ipr", [
    ((), ("IPR011584",)),                       # GFP: a closed 11-strand barrel, soluble
    (("KW-0964",), ("IPR000566",)),             # a secreted lipocalin
    (("KW-0998", "KW-0449"), ("IPR017687",)),   # BamB: a lipoprotein on the outer membrane (beta-propeller)
    (("KW-1134", "KW-0964", "KW-0800"), ()),    # aerolysin: a secreted toxin, its pore forms only in a host
    (("KW-0472", "KW-0812", "KW-1133"), ()),    # a helical membrane protein with a soluble barrel domain
])
def test_soluble_barrels_never_count_as_membrane_barrels(kws, ipr):
    from mpstruc.membrane import counted_barrels, outer_membrane_proteins, raw_facts
    ext, fam, det, cov, feats = _structure("P1", kws, ipr)
    assert outer_membrane_proteins(ext["kw"], fam["ipr"]) == set()
    # Even if a count is handed in, it is not read as a membrane barrel ...
    f = raw_facts("9SOL", ext, fam, det, cov, feats, {"P1": _barrel(11)})
    assert f["barrel_counted"] == {} and not f["barrel_modelled"] and not f["has_barrel"]
    # ... and no coordinates are fetched for such a protein in the first place.
    assert counted_barrels(["9SOL"], {"9SOL": ext}, {"9SOL": fam}, feats, {"9SOL": det}) == {"9SOL": {}}


def test_outer_membrane_evidence_and_only_unmapped_proteins_are_counted(monkeypatch):
    from mpstruc import membrane
    assert membrane.outer_membrane_proteins({"A": ["KW-1134"], "B": ["KW-0998"], "C": ["KW-0998", "KW-0449", "KW-1134"],
                                             "D": [], "E": ["KW-0964", "KW-0998", "KW-1134"]},
                                            {"D": ["IPR039910"], "F": ["IPR000001"]}) == {"A", "B", "C", "D", "E"}
    asked = {}

    def fake(wanted):
        asked.update(wanted)
        return {p: {"chains": {"A": {"entity": "2", "count": 16, "closed": False, "strands": [[404, 410], [780, 786]]}},
                    "oligomers": []} for p in wanted}

    monkeypatch.setattr(membrane.barrel_mod, "barrel_counts", fake)
    ent = lambda *accs: {"polymer_entities": [{"rcsb_id": f"9XYZ_{n}", "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": [a]}}
                                               for n, a in enumerate(accs, 1)]}
    ext = {"9XYZ": {"kw": {"BAMD": ["KW-0998", "KW-0449"], "BAMA": ["KW-0998", "KW-1134"]}}, "9POR": {"kw": {"OMPF": ["KW-0998", "KW-1134"]}}}
    fam = {"9XYZ": {"ipr": {"BAMA": ["IPR023707"]}}, "9POR": {"ipr": {}}}
    seg = {"BAMA": {"tm": [], "named": [[24, 91, "POTRA"]]}, "OMPF": {"tm": [], "named": [[23, 362, "Porin"]]}, "BAMD": {"tm": [], "named": []}}
    ents = {"9XYZ": ent("BAMD", "BAMA"), "9POR": ent("OMPF")}
    cov = {"9XYZ": [{"align": {"BAMD": [[1, 1, 240]]}}, {"align": {"BAMA": [[1, 21, 790]]}}], "9POR": [{"align": {"OMPF": [[1, 23, 340]]}}]}
    got = membrane.counted_barrels(["9XYZ", "9POR"], ext, fam, seg, ents, cov)
    assert asked == {"9XYZ": ["2"]}  # BamA's entity only; the porin has a named barrel domain already
    assert got == {"9XYZ": {"BAMA": {"chain": "A", "count": 16, "method": barrel.METHOD, "closed": False, "entity": "2",
                                     "first": 404, "last": 786}}, "9POR": {}}  # the porin: nothing to count, so no barrel counted
    asked.clear()
    membrane.counted_barrels(["9XYZ", "9POR"], ext, fam, None, ents, cov, seg)  # for the placement rules: every outer-membrane protein
    assert asked == {"9XYZ": ["2"], "9POR": ["1"]}


def test_a_structure_that_could_not_be_counted_is_not_a_structure_without_a_barrel(monkeypatch):
    from mpstruc import membrane
    ent = {"polymer_entities": [{"rcsb_id": "9XYZ_1", "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["BAMA"]}}]}
    ext = {p: {"kw": {"BAMA": ["KW-0998", "KW-1134"]}} for p in ("9DWN", "9BIG", "9NON", "9SOL")}
    ext["9SOL"] = {"kw": {"GFP": []}}
    fam = {p: {"ipr": {}} for p in ext}
    answers = {"9DWN": None,  # the download failed, or the PDB cannot serve it
               "9BIG": {"chains": {}, "oligomers": [], "skipped": "too large"},
               "9NON": {"chains": {}, "oligomers": []}}
    monkeypatch.setattr(membrane.barrel_mod, "barrel_counts", lambda wanted: {p: answers[p] for p in wanted})
    got = membrane.counted_barrels(list(ext), ext, fam, None, {p: ent for p in ext}, {}, {})
    assert got == {"9NON": {}, "9SOL": {}}  # counted, no barrel; 9DWN and 9BIG are absent: not counted


def test_an_open_barrel_with_too_few_strands_moves_nothing():
    # BamA caught with a substrate in it: 13 strands seen, its seam open. 13 would "fit" the
    # autotransporters' 12 +/- 1, but strands may simply be missing from the count.
    bam = _cand("9SOW", "parent", OMP, rule="new-protein", uniprot=("BAMA2",), membrane=0.84)
    twin = _cand("9SOX", "related", OMP, "9SOW", rule="same-paper-batch", uniprot=("BAMA2",))
    # ... while too many strands in an open barrel still contradict a rule (TamA: 16, open).
    tam = _cand("9TAM", "parent", AT, rule="new-protein", uniprot=("TAMA_NEW",))
    barrels = {"9SOW": {"BAMA2": _barrel(13, closed=False)}, "9SOX": {"BAMA2": _barrel(14, closed=False)},
               "9TAM": {"TAMA_NEW": _barrel(16, closed=False)}}
    changes = apply_rules([bam, twin, tam], _kb(), RULED, barrels, {}, {}, {})
    assert (bam["placement"], bam["subgroup"]) == ("parent", OMP) and bam["ruleCheck"]["status"] == "unclear"
    assert "seam open" in bam["ruleCheck"]["text"] and bam["barrelStrands"]["closed"] is False
    assert twin["subgroup"] == OMP and twin["ruleCheck"]["status"] == "unclear"
    assert tam["subgroup"] == OMP and tam["ruleCheck"]["status"] == "moved"
    assert [c["pdb"] for c in changes] == ["9TAM"]
    assert rules.verdict(RULED[OMP], 13, closed=False) == "unclear" and rules.verdict(RULED[OMP], 13, closed=True) == "differs"
    assert rules.verdict(RULED[AT], 16, closed=False) == "differs" and rules.verdict(RULED[OMP], 15, closed=False) == "fits"


def test_a_counted_barrel_never_lowers_a_score():
    from mpstruc.membrane import TM_FLOOR, score_facts
    ext, fam, det, cov, feats = _structure(named=[(24, 91, "POTRA")])
    counted = {"BAMA": _barrel(16)}

    class Ctx:
        tcdb, protein_dates, protein_example, fam_stats = {}, {}, {}, {}
    low = {"bias": -4.0, "weights": {}, "tm_floor": TM_FLOOR, "prior": 0.1}
    f, x, _, p = score_facts("9PXG", ext, fam, det, cov, feats, counted, Ctx, low)
    assert p == TM_FLOOR and f["barrel_counted"] == {"BAMA": 16}
    # A model in which "a barrel is resolved" weighs less than what it replaces: the score stays.
    high = {"bias": 3.0, "weights": {"barrel_modelled": -1.0, "tm_modelled": -1.0}, "tm_floor": TM_FLOOR, "prior": 0.1}
    f, x, _, p = score_facts("9PXG", ext, fam, det, cov, feats, counted, Ctx, high)
    assert f["barrel_counted"] == {} and p == score_facts("9PXG", ext, fam, det, cov, feats, None, Ctx, high)[3] > 0.9
    assert score_facts("9XXX", None, None, None, None, feats, counted, Ctx, low) is None


def test_official_entries_that_contradict_a_rule_are_listed(monkeypatch):
    from types import SimpleNamespace as NS
    from mpstruc import automation, membrane
    entry = lambda pdb, key, name, daughters=(): NS(pdbCode=pdb, key=key, name=name, daughters=list(daughters))
    db = NS(groups=[NS(subgroups=[
        NS(key=AT, parents=[entry("4C00", f"{AT}-P0011", "TamA Autotransporter, full length"), entry("2QOM", f"{AT}-P0003", "EspP"),
                            entry("1WXR", f"{AT}-P0005", "Hbp passenger domain")]),
        NS(key=OMP, parents=[entry("5D0O", f"{OMP}-P0010", "BamABCDE", [entry("7TTC", f"{OMP}-P0010-D001", "BamABCDE with EspP")]),
                             entry("9HJM", f"{OMP}-P0021", "BamADGHIJK")]),
        NS(key="G03-S001", parents=[entry("1ABC", "G03-S001-P0001", "a channel")])])])
    for name in ("entry_extras", "entry_families", "entry_entities"):
        monkeypatch.setattr(automation.rcsb, name, lambda ids: {p: {} for p in ids})
    found = {"4C00": {"TAMA": _barrel(16)}, "2QOM": {"ESPP": _barrel(12)}, "5D0O": {"BAMA": _barrel(16, closed=False)},
             "7TTC": {"ESPP": _barrel(12, "F"), "BAMA": _barrel(15, closed=False)}, "9HJM": {"BAMA": _barrel(14, closed=False)}}
    asked = []
    monkeypatch.setattr(membrane, "counted_barrels", lambda ids, *a: asked.extend(ids) or found)
    got = automation.rule_conflicts(db, _kb(), RULED)
    assert sorted(asked) == ["1WXR", "2QOM", "4C00", "5D0O", "7TTC", "9HJM"]  # only the ruled subfamilies
    assert [(c["pdb"], c["count"], c["expected"], c["fits"]) for c in got] == [("4C00", 16, 12, OMP)]
    md = automation.conflicts_md(got, RULED)
    assert "| 4C00 | TamA Autotransporter, full length" in md and "Typically 12-stranded beta-barrel." in md
    assert "None." in automation.conflicts_md([], RULED)


# ---------------------------------------------------------------- review findings

def test_an_open_count_is_a_lower_bound_for_every_rule():
    at, omp = RULED[AT], RULED[OMP]
    # 11 to 13 strands in an open barrel are in the autotransporters' range, but a 16-strand
    # Omp85 barrel caught open looks the same: the count decides nothing.
    for n in (11, 12, 13):
        assert rules.verdict(at, n, closed=False, others=[omp]) == "unclear"
        assert rules.verdict(at, n, closed=True, others=[omp]) == "fits"
    assert rules.verdict(at, 14, closed=False, others=[omp]) == "differs"      # more than 12 +/- 1 whatever is missing
    assert rules.verdict(omp, 16, closed=False, others=[at]) == "fits"         # no rule with more strands is in reach
    assert rules.verdict(omp, 9, closed=False, others=[at]) == "unclear"
    assert rules.possible(omp, 13, closed=False) and not rules.possible(at, 14, closed=False)
    assert rules.possible(omp, 15) and not rules.possible(omp, 13)
    # Two rings fit the same strands (8, or 10 with a hairpin): undecided if only one number fits.
    ten = {"barrelStrands": 10, "tolerance": 0}
    assert rules.verdict(ten, 8, up_to=10) == "unclear" and rules.verdict(at, 8, up_to=10) == "differs"
    assert barrel.strands_text(_barrel(16)) == "16" and barrel.strands_text(_barrel(13, closed=False)) == "at least 13"
    assert barrel.strands_text({**_barrel(8), "upTo": 10}) == "8 to 10"
    assert barrel.barrel_text(_barrel(8)) == "an 8-strand beta-barrel"
    assert barrel.barrel_text(_barrel(13, closed=False)) == "a beta-barrel of at least 13 strands"


def test_an_open_undercount_is_not_confirmed_as_an_autotransporter():
    # A new TamA/BamA-like structure that family votes put under Autotransporters, counted 13 open.
    c = _cand("9OPN", "parent", AT, rule="new-protein", uniprot=("OMP85_NEW",))
    assert apply_rules([c], _kb(), RULED, {"9OPN": {"OMP85_NEW": _barrel(13, closed=False)}}, {"9OPN": {"OMP85_NEW"}}, {}, {}) == []
    assert c["ruleCheck"]["status"] == "unclear" and c["subgroup"] == AT and c["placement"] == "parent"
    assert c["ruleCheck"]["text"].startswith("At least 13 strands counted in the barrel (chain A, seam open): that would fit")
    assert "Omp85-TpsB Outer Membrane Transporter Superfamily are typically 16-stranded" in c["ruleCheck"]["text"]
    # 14 open: too many for an autotransporter, not yet shown to be 16: noted, not moved.
    d = _cand("9OPM", "parent", AT, rule="new-protein", uniprot=("OMP85_NEW",))
    apply_rules([d], _kb(), RULED, {"9OPM": {"OMP85_NEW": _barrel(14, closed=False)}}, {"9OPM": {"OMP85_NEW"}}, {}, {})
    assert d["ruleCheck"]["status"] == "differs" and d["subgroup"] == AT and "it may be one of Omp85-TpsB" in d["ruleCheck"]["text"]


def _bam(pdb="9PXG", membrane=0.45):
    return _cand(pdb, "parent", OMP, rule="new-protein", uniprot=("BAMA2",), membrane=membrane)


def test_an_uncounted_bam_complex_is_not_sent_to_review():
    # The coordinate server is down (rule_inputs returns no barrels at all), or this one
    # download failed: the structure has an outer-membrane protein nobody looked at.
    for barrels in ({}, {"9OTH": {}}):
        lead, twin = _bam(), _cand("9PXJ", "related", None, "9PXG", rule="same-paper-batch", uniprot=("BAMA2",), membrane=0.48)
        before = json.loads(json.dumps([lead, twin]))
        assert apply_rules([lead, twin], _kb(), RULED, barrels, {"9PXG": {"BAMA2"}, "9PXJ": {"BAMA2"}}, {}, {}) == []
        assert lead["ruleCheck"]["status"] == "unchecked" and "could not be counted" in lead["ruleCheck"]["text"]
        for now, was in zip((lead, twin), before):
            assert {k: now.get(k) for k in ("placement", "subgroup", "ownerPdb", "rule", "confidence", "evidence")} == \
                   {k: was.get(k) for k in ("placement", "subgroup", "ownerPdb", "rule", "confidence", "evidence")}
        assert "ruleCheck" not in twin and "barrelStrands" not in lead
    # Counted, and no barrel in it: that is what sends it to review.
    lead = _bam()
    apply_rules([lead], _kb(), RULED, {"9PXG": {}}, {"9PXG": {"BAMA2"}}, {}, {})
    assert lead["placement"] == "review" and lead["ruleCheck"]["status"] == "no-barrel"
    # A structure with no outer-membrane protein has nothing to count: the pore still goes to review.
    pore = _cand("9AAV", "parent", AT, rule="new-protein", uniprot=("NUP",), membrane=0.171)
    apply_rules([pore], _kb(), RULED, {}, {}, {}, {})
    assert pore["placement"] == "review" and pore["subgroup"] is None


CONFLICT_4C00 = [{"pdb": "4C00", "key": f"{AT}-P0011", "subgroup": AT, "count": 16, "chain": "A", "closed": True,
                  "uniprot": "TAMA_ECOLI", "fits": OMP}]


def _tam_paper():
    return [_cand("9TGQ", "review", AT, "4C00"), _cand("9TGR", "related", AT, "9TGQ", rule="same-paper-batch"),
            _cand("9TGS", "related", AT, "9TGQ", rule="same-paper-batch")]


def test_the_tama_paper_moves_even_if_the_lead_has_no_usable_barrel():
    fams = {"9TGQ": {"ipr": {"TAMA_ECOLI": ["IPR000184", "IPR035243"]}}}
    om = {p: {"TAMA_ECOLI"} for p in ("9TGQ", "9TGR", "9TGS")}
    # 1. No barrel found in the lead (6 A map), 16 closed strands in a structure of the same paper.
    props = _tam_paper()
    apply_rules(props, _kb(), RULED, {"9TGQ": {}, "9TGR": {}, "9TGS": {"TAMA_ECOLI": _barrel(16)}}, om, fams, {})
    lead, r, s = props
    assert lead["subgroup"] == OMP and lead["ruleCheck"]["status"] == "moved" and lead["ownerPdb"] == "8US1"
    assert lead["ruleCheck"]["text"].startswith("16 strands counted in the barrel (chain A of 9TGS, same paper and protein): Omp85")
    assert "barrelStrands" not in lead  # nothing was counted in the lead itself
    assert r["subgroup"] == s["subgroup"] == OMP and s["barrelStrands"]["count"] == 16 and "barrelStrands" not in r
    # 2. The lead's own barrel is open and short (14: more than 12, not shown to be 16); same follower.
    props = _tam_paper()
    apply_rules(props, _kb(), RULED, {"9TGQ": {"TAMA_ECOLI": _barrel(14, closed=False)}, "9TGS": {"TAMA_ECOLI": _barrel(16)}}, om, fams, {})
    assert props[0]["subgroup"] == OMP and props[0]["barrelStrands"]["count"] == 14 and "9TGS" in props[0]["ruleCheck"]["text"]
    # 3. Nothing counted anywhere, but the official entry it would join (4C00) is itself against the rule.
    props = _tam_paper()
    apply_rules(props, _kb(), RULED, {"9TGQ": {}, "9TGR": {}, "9TGS": {}}, om, fams, {}, CONFLICT_4C00)
    lead = props[0]
    assert lead["subgroup"] == OMP and lead["ruleCheck"]["status"] == "moved" and lead["ruleCheck"]["from"] == AT
    assert "16 strands are counted in the barrel of 4C00 itself" in lead["ruleCheck"]["text"]
    assert all(c["subgroup"] == OMP for c in props)
    # ... also when the coordinates could not be fetched at all (not counted).
    props = _tam_paper()
    apply_rules(props, _kb(), RULED, {}, om, fams, {}, CONFLICT_4C00)
    assert props[0]["subgroup"] == OMP and props[0]["ruleCheck"]["status"] == "moved"
    # 4. No count, no conflict list: UniProt's 16 strands, all resolved (the per-protein table), decide.
    props = _tam_paper()
    props[0]["tmEvidence"] = [{"uniprot": "TAMA_ECOLI", "segments": 16, "resolved": 16, "status": "resolved"},
                              {"uniprot": "TAMB", "segments": 1, "resolved": 0, "status": "unresolved"}]
    apply_rules(props, _kb(), RULED, {"9TGQ": {}, "9TGR": {}, "9TGS": {}}, om, fams, {})
    lead = props[0]
    assert lead["subgroup"] == OMP and lead["barrelStrands"] == {"count": 16, "method": "UniProt strands", "closed": True, "uniprot": "TAMA_ECOLI"}
    assert lead["ruleCheck"]["text"].startswith("All 16 membrane-crossing strands UniProt gives for TAMA_ECOLI are resolved")
    # 5. The owner's conflict never overrides a barrel of the structure that fits where it is.
    twelve = _cand("9FIT", "review", AT, "4C00")
    apply_rules([twelve], _kb(), RULED, {"9FIT": {"TAMA_ECOLI": _barrel(12)}}, {"9FIT": {"TAMA_ECOLI"}}, fams, {}, CONFLICT_4C00)
    assert twelve["subgroup"] == AT and twelve["ruleCheck"]["status"] == "fits"


def test_a_ring_made_by_several_chains_never_moves_a_suggestion_to_omp85():
    # (What a head domain of a trimeric autotransporter looked like before the anchor test;
    # Omp85 and TpsB barrels are made by one chain.)
    c = _cand("9TAA", "parent", AT, rule="new-protein", uniprot=("BADA",))
    b = {**_barrel(15), "chains": 3, "perChain": 5.0}
    assert apply_rules([c], _kb(), RULED, {"9TAA": {"BADA": b}}, {"9TAA": {"BADA"}}, {}, {}) == []
    assert (c["subgroup"], c["placement"], c["ruleCheck"]["status"]) == (AT, "parent", "differs")
    assert "made by one chain, so it is not moved" in c["ruleCheck"]["text"]


def test_a_barrel_in_a_protein_uniprot_does_not_place_in_the_outer_membrane_only_confirms():
    # Invasin: no outer-membrane keyword, so its chains are counted only because the
    # structure is suggested under a ruled subfamily ("vote"). 12 closed strands confirm it ...
    inv = _cand("9INV", "parent", AT, rule="new-protein", uniprot=("INVASIN",), membrane=0.3)
    apply_rules([inv], _kb(), RULED, {"9INV": {"INVASIN": {**_barrel(12), "vote": True}}}, {}, {}, {})
    assert (inv["placement"], inv["subgroup"], inv["ruleCheck"]["status"]) == ("parent", AT, "fits")
    # ... but a soluble barrel of another size (a lipocalin fused to something) is not evidence of anything.
    sol = _cand("9LIP", "parent", AT, rule="new-protein", uniprot=("LIPOCALIN",), membrane=0.3)
    apply_rules([sol], _kb(), RULED, {"9LIP": {"LIPOCALIN": {**_barrel(8), "vote": True}}}, {}, {}, {})
    assert (sol["placement"], sol["subgroup"], sol["ruleCheck"]["status"]) == ("review", None, "no-barrel")


# The membrane parts UniProt names for a trimeric autotransporter (Hia: only a named
# anchor region; BadA: strands and a translocator domain) and where a construct's ring is.
HIA = {"tm": [], "named": [[53, 104, "Trimeric autotransporter adhesin putative GIN"],
                           [1038, 1098, "Trimeric autotransporter adhesin YadA-like C-terminal membrane anchor"]]}
BADA = {"tm": [[3028, 3039], [3044, 3051], [3056, 3066], [3070, 3082]],
        "named": [[48, 2901, "Surface exposed passenger domain"], [2902, 3027, "Outer membrane translocation of the passenger domain"],
                  [3028, 3082, "Translocator domain"]]}


def _trimer(count, first, last):
    return {"chain": "A", "count": count, "method": barrel.METHOD, "closed": True, "entity": "1", "chains": 3,
            "perChain": round(count / 3, 1), "first": first, "last": last}


@pytest.mark.parametrize("name, feat, align, b, want", [
    ("Hia anchor 2GR8", HIA, [[23, 1022, 77]], _trimer(12, 40, 95), True),
    ("Hia binding domain 1S7M", HIA, [[1, 541, 174]], _trimer(8, 113, 165), False),
    ("BadA head 3D9X: 15 strands, 5 per chain", BADA, [[1, 374, 163]], _trimer(15, 20, 150), False),
    ("AtaA-like stalk fragment, 17 strands", BADA, [[30, 2000, 163]], _trimer(17, 40, 180), False),
    ("a translocator domain alone", BADA, [[1, 3000, 83]], _trimer(12, 30, 80), True),
    ("the chain cannot be put in UniProt numbering", BADA, [], _trimer(12, 30, 80), False),
    ("MspA: UniProt names no membrane part, 8 x 2 strands", {"tm": [], "named": []}, [[1, 28, 184]],
     {**_trimer(16, 36, 183), "chains": 8, "perChain": 2.0}, True),
    ("no membrane part named, 8 strands over 3 chains", {"tm": [], "named": []}, [[1, 28, 184]], _trimer(8, 36, 183), False),
    ("one chain, no membrane part named (BamA)", {"tm": [], "named": [[24, 91, "POTRA"]]}, [[1, 21, 790]],
     {**_barrel(16), "first": 404, "last": 786}, True),
    ("one chain, ring outside the strands UniProt gives", BADA, [[1, 374, 400]], {**_barrel(9), "first": 20, "last": 150}, False),
])
def test_only_a_ring_in_the_membrane_part_is_a_membrane_barrel(name, feat, align, b, want):
    from mpstruc.membrane import membrane_barrel
    cov = [{"align": {"P1": align}, "unobs": []}]
    assert membrane_barrel(b, "P1", cov, feat, ("KW-0998", "KW-1134")) is want, name


def test_soluble_heads_of_trimeric_autotransporters_do_not_raise_the_score_or_touch_the_rules(monkeypatch):
    from mpstruc import membrane
    ring = {"chains": {}, "oligomers": [{"entity": "1", "chains": ["A", "B", "C"], "count": 15, "perChain": 5.0, "first": 20, "last": 150}]}
    monkeypatch.setattr(membrane.barrel_mod, "barrel_counts", lambda wanted: {p: ring for p in wanted})
    ent = {"polymer_entities": [{"rcsb_id": "3D9X_1", "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["BADA"]}}]}
    ext, fam = {"3D9X": {"kw": {"BADA": ["KW-0998", "KW-1134"]}}}, {"3D9X": {"ipr": {}}}
    head = {"3D9X": [{"align": {"BADA": [[1, 374, 163]]}, "unobs": []}]}
    anchor = {"3D9X": [{"align": {"BADA": [[1, 2940, 143]]}, "unobs": []}]}
    assert membrane.counted_barrels(["3D9X"], ext, fam, None, {"3D9X": ent}, head, {"BADA": BADA}) == {"3D9X": {}}
    got = membrane.counted_barrels(["3D9X"], ext, fam, None, {"3D9X": ent}, anchor, {"BADA": BADA})
    assert got["3D9X"]["BADA"]["count"] == 15 and got["3D9X"]["BADA"]["chains"] == 3


def test_the_guard_set_in_the_validation_report_is_clean():
    # reports/barrel_validation.json is written by `python -m mpstruc.barrel validate` from
    # the real structures: no soluble barrel or soluble fragment of an outer-membrane
    # protein may be read as a membrane barrel, and the report must cover the whole set.
    from mpstruc.report import ROOT
    v = json.loads((ROOT / "reports" / "barrel_validation.json").read_text())
    assert v["guard_failed"] == []
    assert {r["pdb"] for r in v["guard"]} == {p for p, _ in barrel.GUARD} >= {"1S7M", "3LA9", "4USX", "3D9X", "3WPR"}
    assert {r["pdb"] for r in v["soluble"]} == {p for p, _ in barrel.SOLUBLE}
    assert all(r["membrane"] == "no" for r in v["guard"] + v["soluble"])
    assert v["open_n"] >= 10 and {r["pdb"] for r in v["open"]} == {p for p, *_ in barrel.OPEN_SEAM}


def test_the_assembly_is_searched_even_when_the_deposited_chains_make_a_small_ring():
    # MspA 1UUN: two of eight chains deposited; a ring of 16 only in the assembly.
    rot = lambda a: [[[math.cos(a), -math.sin(a), 0.0], [math.sin(a), math.cos(a), 0.0], [0.0, 0.0, 1.0]], [0.0, 0.0, 0.0]]
    chains = {c: _ring(2, slots=16, first=0, phase=2 * math.pi * 2 * k / 16) for k, c in enumerate("AB")}
    ops = [rot(2 * math.pi * k / 4) for k in range(4)]
    res = barrel.count_barrels({"chains": chains, "operators": ops})
    (o,) = res["oligomers"]
    assert (o["count"], len(o["chains"]), o["perChain"]) == (16, 8, 2.0) and (o["first"], o["last"]) == (1, 24)
    # Without the operators the two chains make no barrel.
    assert barrel.count_barrels({"chains": chains, "operators": []})["oligomers"] == []


def test_a_ring_that_two_more_strands_could_join_is_reported_as_a_range():
    # Chlamydia MOMP 9OQ7: the last and the first strand of an 8-ring are bridged (9 pairs),
    # and a hairpin after the last strand bridges both more strongly (20 and 18 pairs).
    ring = list(range(8))
    edges = {(k, k + 1): 15 for k in range(7)}
    edges[(0, 7)] = 9
    assert barrel.ring_detours(ring, set(range(8)), edges) == 0
    edges.update({(7, 8): 20, (8, 9): 21, (0, 9): 18})
    assert barrel.ring_detours(ring, set(range(10)), edges) == 2
    edges[(0, 9)] = 5  # a hairpin that only leans on the last strand: not a second ring
    assert barrel.ring_detours(ring, set(range(10)), edges) == 0
    b = barrel.best_barrel({"chains": {"A": {"entity": "1", "count": 8, "closed": True, "upTo": 10, "strands": [[1, 9], [90, 99]]}}, "oligomers": []})
    assert b["upTo"] == 10 and barrel.strands_text(b) == "8 to 10"


def test_subfamily_votes_name_the_same_features_every_run():
    from mpstruc.automation import subgroup_scores
    kb = _kb()
    feats = [f"ipr:IPR{n:06d}" for n in range(1, 9)] + ["uipr:IPR000184", "pfam:PF01103"]
    kb.feat_sg = {f: {OMP: 3, AT: 1} for f in feats}
    s = {"uniprot": set(), "pfam": {"PF01103"}, "interpro": {f.split(":")[1] for f in feats[:8]}, "opm": set(), "uipr": {"IPR000184"}, "core": set()}
    runs = [subgroup_scores(kb, json.loads(json.dumps({k: sorted(v, reverse=bool(n % 2)) for k, v in s.items()}))) for n in range(2)]
    assert runs[0] == runs[1] and [sg for sg, _, _ in runs[0]] == [OMP, AT]
    assert len(runs[0][0][2]) == 4 and runs[0][0][2][0].split(":")[0] in ("pfam", "uipr")  # the heaviest features first
