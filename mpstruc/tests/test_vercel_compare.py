import csv
import json

from mpstruc import vercel_compare as vc
from mpstruc.vercel_compare import Sources, compare_week, compare_whole, judge, merge, release_wednesday


def _snap():
    e = lambda status, fam=None: {"status": status, "family": fam, "name": None, "group": None, "released": None,
                                  "tm": None, "memo": None, "seenAt": "2026-11-13"}
    return {
        "1AAA": e("already_in", "Multi-Pass"),  # in mpstruc
        "2BBB": e("no", "non-membrane"),  # our model says not membrane; Vercel said No too
        "3CCC": e("pending", "Peripheral"),  # our curators rejected it as not membrane
        "4DDD": e("already_in", "Multi-Pass"),  # 'Already in' but not in mpstruc
        "5EEE": e("yes", "Single-Pass"),  # model uncertain
    }


def _src(**kw):
    scores = {"2BBB": {"p": 0.02, "why": ["UniProt does not list any of its proteins as crossing the membrane"],
                       "facts": {"has_tm_protein": False, "has_barrel": False, "tm_in_construct": 0, "tm_modelled": 0}},
              "5EEE": {"p": 0.3, "why": ["the title mentions 'receptor'"], "facts": None},
              "6FFF": {"p": 0.92, "why": ["11 membrane-crossing segments are resolved in the structure",
                                          "the sample contains detergent (LMN)"],
                       "facts": {"has_tm_protein": True, "tm_in_construct": 12, "tm_modelled": 11}},
              "7GGG": {"p": 0.6, "why": ["UniProt lists a protein in it as crossing the membrane"], "facts": None},
              "9III": {"p": 0.05, "why": [], "facts": {"has_tm_protein": True, "tm_in_construct": 4, "tm_modelled": 0}}}
    asked = []

    def scorer(ids):
        asked.extend(ids)
        return {p: scores.get(p) for p in ids}
    src = Sources(
        curation={"3CCC": {"status": "rejected", "reasonCode": "not_membrane", "reason": "Not a membrane protein"},
                  "8HHH": {"status": "rejected", "reasonCode": "unpublished", "reason": "Not published yet"},
                  "1ZZZ": {"status": "approved"}},
        mpstruc={"1AAA": "Listed in the official mpstruc (Channels: Potassium)",
                 "3CCC": "Listed in the official mpstruc (Monotopic)",
                 "8HHH": "Listed in the official mpstruc as related to 1AAA (Channels: Potassium)",
                 "1YYY": "Listed in the official mpstruc (Porins)"},
        mpstruc_entries={"1AAA", "3CCC", "1YYY"},
        pdbtm={"8HHH"}, opm={"7GGG": "monotopic/peripheral", "1XXX": "transmembrane"},
        ours={"6FFF", "3CCC"}, scorer=scorer, **kw)
    return src, asked


def test_evidence_precedence():
    src, _ = _src()
    src.score(["2BBB", "5EEE", "6FFF", "7GGG", "9III"])
    # Curators' 'not a membrane protein' beats the mpstruc listing.
    v, ev, _, by = judge("3CCC", src)
    assert v == "not_membrane" and by == "curators" and ev[0]["text"] == "Our curators rejected it: Not a membrane protein"
    assert ev[1]["source"] == "mpstruc"
    # A rejection for another reason does not decide; mpstruc does.
    v, ev, _, by = judge("8HHH", src)
    assert v == "membrane" and by == "mpstruc" and "says nothing about the membrane" in ev[0]["text"]
    assert any(e["source"] == "pdbtm" for e in ev)
    assert judge("1ZZZ", src)[0] == "membrane"
    assert judge("1XXX", src)[:2] == ("membrane", [{"source": "opm", "text": vc.OPM_TEXT["transmembrane"]}])
    # OPM monotopic alone does not decide; the model does.
    v, ev, _, by = judge("7GGG", src)
    assert v == "membrane" and by == "model" and [e["source"] for e in ev] == ["opm", "model"]
    v, ev, _, _ = judge("6FFF", src)
    assert v == "membrane" and ev[0]["text"].startswith("Our model: 92% - 11 membrane-crossing segments resolved (of 12")
    assert ev[0]["text"].count("membrane-crossing segment") == 1  # the model's own copy of the fact is not repeated
    v, ev, _, _ = judge("2BBB", src)
    assert v == "not_membrane" and "UniProt marks no membrane-crossing segment" in ev[0]["text"]
    assert "none of the 4 membrane-crossing segments in the construct is resolved" in judge("9III", src)[1][0]["text"]
    assert judge("5EEE", src)[0] == "uncertain"
    assert judge("0NON", src) == ("uncertain", [], 0.0, "none")


def test_compare_week_missed_and_not_membrane():
    src, asked = _src()
    snap = _snap()
    ids = ["1AAA", "2BBB", "3CCC", "5EEE", "6FFF", "7GGG", "9III", "1XXX"]
    w = compare_week("2026-11-11", ids, src, snap)
    assert w["week"] == "2026-11-11" and w["from"] == "2026-11-05"
    # Curators, PDBTM and OPM (with segments) need no score; mpstruc listings are checked against the facts.
    assert set(asked) == {"1AAA", "2BBB", "5EEE", "6FFF", "7GGG", "9III"}
    assert w["vercelListed"] == 4  # 1AAA 2BBB 3CCC 5EEE
    assert w["agree"] == 1  # 1AAA
    assert [r["pdb"] for r in w["missed"]] == ["1XXX", "6FFF", "7GGG"]  # databases first, then by score
    assert [r["pdb"] for r in w["notMembrane"]] == ["3CCC", "2BBB"]  # curators first
    r = w["notMembrane"][1]
    assert r["vercelStatus"] == "no" and r["vercelFamily"] == "non-membrane" and r["ourScore"] == 0.02
    assert w["missed"][1]["inOurList"] and not w["missed"][1]["inMpstruc"]
    assert set(w["missed"][0]) == {"pdb", "title", "released", "verdict", "evidence", "ourScore", "vercelStatus",
                                   "vercelFamily", "inMpstruc", "inOurList"}
    assert w["verdicts"] == {"membrane": 4, "not_membrane": 3, "uncertain": 1}


def test_compare_week_counts_page_read_entries_by_release_date():
    src, _ = _src()
    snap = {"0PGE": {"status": "pending", "released": "2026-11-11"}}
    assert compare_week("2026-11-11", ["6FFF"], src, snap)["vercelListed"] == 1


def test_compare_whole():
    src, asked = _src()
    snap = _snap()
    span = ["1AAA", "1YYY", "1XXX", "6FFF", "7GGG", "2BBB"]
    whole = compare_whole(snap, src, span, recent_ids=["6FFF", "7GGG"])
    s = whole["summary"]
    assert [r["pdb"] for r in whole["notMembrane"]] == ["3CCC", "2BBB"]
    assert s["notMembraneVercelSaidNo"] == 1 and s["notMembraneByVercelStatus"] == {"pending": 1, "no": 1}
    # Strong sources over the span; model only for recent releases.
    assert [r["pdb"] for r in whole["missed"]] == ["1YYY", "1XXX", "6FFF", "7GGG"]
    assert s["missedStrongSources"] == 2 and s["missedModelRecent"] == 2
    assert sum(s["missedBySource"].values()) == s["missed"] and s["missedBySource"].get("model") == 2
    assert s["alreadyInNotInMpstruc"] == {"count": 1, "pdbs": ["4DDD"]}
    m = s["mpstrucNotInVercel"]
    assert m["pdbs"] == ["1YYY", "8HHH"] and m["entries"] == 1 and m["relatedOnly"] == 1 and m["inVercelSpan"] == 1
    assert s["decidedBy"] == {"mpstruc": 1, "model": 1, "curators": 1, "none": 2}
    assert "3CCC" not in asked and "1AAA" in asked  # curators decide 3CCC; the mpstruc listing is checked


def test_merge_keeps_past_weeks_and_whole():
    old = {"weeks": [{"week": "2026-11-04", "agree": 1}, {"week": "2026-11-11", "agree": 0}],
           "wholeDatabase": {"summary": {"x": 1}, "missed": [], "notMembrane": []}}
    out = merge(old, {"entries": 3}, week={"week": "2026-11-11", "agree": 5})
    assert [(w["week"], w["agree"]) for w in out["weeks"]] == [("2026-11-11", 5), ("2026-11-04", 1)]
    assert out["wholeDatabase"]["summary"] == {"x": 1} and out["snapshot"] == {"entries": 3}
    out = merge(out, {}, whole={"summary": {"y": 2}, "missed": [], "notMembrane": []})
    assert out["wholeDatabase"]["summary"] == {"y": 2} and len(out["weeks"]) == 2


def test_release_wednesday():
    import datetime as dt
    assert release_wednesday(dt.date(2026, 11, 13)).isoformat() == "2026-11-11"
    assert release_wednesday(dt.date(2026, 11, 11)).isoformat() == "2026-11-11"
    assert release_wednesday(dt.date(2026, 11, 10)).isoformat() == "2026-11-04"


def test_snapshot_from_export(tmp_path):
    p = tmp_path / "mpstruc_form_export_2026-11-13.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Pdb Code", "Name", "Family", "Bibliography", "PubMed Link", "Status", "Memo"])
        w.writerow(["9tgs", "TamA", "Outer-Membrane", "Paper. Cell Rep (2026). PMID:42507553", "", "Already in", ""])
        w.writerow(["8ABC", "Thing", "", "", "https://www.ncbi.nlm.nih.gov/pubmed/?term=123", "", "check"])
    s = vc.snapshot_from_export(p)
    assert s["source"] == "export:mpstruc_form_export_2026-11-13.csv"
    assert s["entries"]["9TGS"]["status"] == "already_in" and s["entries"]["9TGS"]["group"] == "42507553"
    assert s["entries"]["8ABC"] == {"status": "pending", "family": None, "name": "Thing", "group": "123",
                                    "released": None, "tm": None, "memo": "check", "seenAt": "2026-11-13"}
    assert vc.snapshot_counts(s)["byStatus"] == {"already_in": 1, "pending": 1}


def test_main_week_end_to_end(tmp_path, monkeypatch):
    snap = tmp_path / "snap.json"
    snap.write_text(json.dumps({"source": "page", "readAt": "2026-11-13T00:00:00Z", "entries": _snap()}))
    out = tmp_path / "out.json"
    out.write_text(json.dumps({"weeks": [{"week": "2026-11-04", "missed": [], "notMembrane": [], "agree": 0,
                                          "vercelListed": 0}]}))
    src, _ = _src(summaries={"6FFF": {"title": "A transporter", "released": "2026-11-11"}})
    monkeypatch.setattr(vc, "load_sources", lambda: src)
    monkeypatch.setattr(vc, "annotated_cached", lambda *a: set())
    import mpstruc.automation as auto
    import mpstruc.rcsb as rcsb
    monkeypatch.setattr(auto, "released_between", lambda a, b: ["1AAA", "2BBB", "6FFF"])
    monkeypatch.setattr(rcsb, "entry_summaries", lambda ids: {})  # no network in tests
    vc.main(["--week", "2026-11-11", "--snapshot", str(snap), "--out", str(out)])
    d = json.loads(out.read_text())
    assert [w["week"] for w in d["weeks"]] == ["2026-11-11", "2026-11-04"]
    assert d["weeks"][0]["missed"][0]["title"] == "A transporter"
    assert d["snapshot"]["entries"] == 5 and d["snapshot"]["source"] == "page"


# ---------------------------------------------------------------- review regressions

def test_opm_needs_a_segment_in_the_structure():
    from mpstruc import membranedb
    detail = {"2JT5": {"type": "transmembrane", "segments": 0, "classtype": 11},  # MMP-3 catalytic domain
              "3V66": {"type": "transmembrane", "segments": 0, "classtype": 11},  # squalene synthase, soluble
              "1D3G": {"type": "transmembrane", "segments": 0, "classtype": 11},  # DHODH
              "2RH1": {"type": "transmembrane", "segments": 14, "classtype": 1},  # beta2 adrenergic receptor
              "1AAP": {"type": "monotopic/peripheral", "segments": 0, "classtype": 2}}
    lab = membranedb.opm_spanning(detail)
    assert lab["2RH1"] == "transmembrane" and lab["1AAP"] == "monotopic/peripheral"
    assert {lab[p] for p in ("2JT5", "3V66", "1D3G")} == {membranedb.OPM_NO_SEGMENT}
    src = Sources(opm=lab, scores={
        "2JT5": {"p": 0.001, "why": [], "facts": {"has_tm_protein": False, "has_barrel": False,
                                                  "tm_in_construct": 0, "tm_modelled": 0, "barrel_modelled": 0}},
        "1D3G": {"p": 0.009, "why": [], "facts": {"has_tm_protein": True, "has_barrel": False,
                                                  "tm_in_construct": 0, "tm_modelled": 0, "barrel_modelled": 0}}})
    assert judge("2RH1", src)[0] == "membrane" and vc.decided_without_model("2RH1", src)
    assert not vc.decided_without_model("2JT5", src)
    v, ev, _, by = judge("2JT5", src)
    assert v == "not_membrane" and by == "model" and "no vote" in ev[0]["text"]
    assert judge("1D3G", src)[0] == "not_membrane"  # DHODH: same verdict as 9OC9 (decided by the model)
    assert judge("3V66", src)[0] == "uncertain"  # no score: OPM alone no longer decides


# Facts our scorer reported for the structures the review found wrongly 'missed'.
_NONE = {"has_tm_protein": False, "has_barrel": False, "tm_in_construct": 0, "tm_modelled": 0,
         "barrel_in_construct": 0, "barrel_modelled": 0}
_CUT = {**_NONE, "has_tm_protein": True}
_REVIEWED = {  # pdb: (p, facts, mpstruc?, opm label)
    "2JT5": (0.001, _NONE, False, "transmembrane-family, no segment in structure"),
    "3V66": (0.185, {**_CUT, "tm_in_construct": 1, "tm_modelled": 1}, False, "transmembrane-family, no segment in structure"),
    "1DYL": (0.007, _CUT, False, "transmembrane-family, no segment in structure"),
    "3NWD": (0.005, _CUT, False, "transmembrane-family, no segment in structure"),
    "8JLE": (0.3, _CUT, True, None),  # SV2A luminal domain + BoNT Hc
    "6IJ6": (0.0, _NONE, True, None),  # PETase at RCSB; mpstruc meant 6I6J
    "6KFI": (0.005, _NONE, True, None),  # G-quadruplex DNA; mpstruc meant 6FKI
    "7JOI": (0.002, _NONE, True, None),  # DNA lattice; mpstruc meant 7OJI
    "7FR3": (0.001, _CUT, True, None),  # SARS-CoV-2 Nsp3 domain; mpstruc meant 7RF3
    "7TM9": (0.001, _NONE, True, None),  # alkaline phosphatase; mpstruc meant 7T9M
}


def _reviewed_src():
    return Sources(
        mpstruc={p: "Listed in the official mpstruc (X)" for p, t in _REVIEWED.items() if t[2]},
        mpstruc_entries={p for p, t in _REVIEWED.items() if t[2]},
        opm={p: t[3] for p, t in _REVIEWED.items() if t[3]},
        scores={p: {"p": t[0], "why": [], "facts": t[1]} for p, t in _REVIEWED.items()})


def test_reviewed_wrong_codes_are_not_missed():
    src = _reviewed_src()
    snap = {"7OJI": {"status": "already_in"}, "7RF3": {"status": "already_in"}, "1AAA": {"status": "no"}}
    whole = compare_whole(snap, src, span_ids=list(_REVIEWED), recent_ids=list(_REVIEWED))
    assert whole["missed"] == []
    w = compare_week("2026-11-11", list(_REVIEWED), src, snap)
    assert w["missed"] == []
    v, ev, _, _ = judge("6IJ6", src)
    assert v == "not_membrane" and "possible mpstruc ID typo" in ev[0]["text"]
    v, ev, _, _ = judge("8JLE", src)
    assert v == "uncertain" and "fragment without the membrane part" in ev[0]["text"]
    s = whole["summary"]
    # 7JOI/7OJI and 7FR3/7RF3 are paired, and left out of both summary lists.
    assert [(x["mpstruc"], x["vercel"]) for x in s["possibleMpstrucTypos"]["pairs"]] == [("7FR3", "7RF3"), ("7JOI", "7OJI")]
    assert s["alreadyInNotInMpstruc"] == {"count": 0, "pdbs": []}
    assert s["mpstrucNotInVercel"]["pdbs"] == []
    w = s["mpstrucWithoutMembranePart"]
    assert w["noMembraneProtein"] == ["6IJ6", "6KFI", "7TM9"] and w["fragments"] == ["8JLE"]


def test_mpstruc_kept_when_membrane_part_resolved_or_unknown():
    src = Sources(mpstruc={p: "Listed in the official mpstruc (X)" for p in ("1OKA", "2UNM", "3NOS")}, scores={
        "1OKA": {"p": 0.2, "why": [], "facts": {**_CUT, "tm_in_construct": 6, "tm_modelled": 6}},
        "2UNM": {"p": 0.05, "why": [], "facts": {**_CUT, "tm_unmapped": ["P0A910"]}}})  # unmapped barrel
    assert [judge(p, src)[3] for p in ("1OKA", "2UNM", "3NOS")] == ["mpstruc"] * 3


def test_memprotmd_alone_needs_some_support():
    src = Sources(memprotmd={"7O7P", "9DIS", "5ABC"}, scores={
        "7O7P": {"p": 0.031, "why": [], "facts": _NONE},  # alpha-2-macroglobulin
        "9DIS": {"p": 0.917, "why": [], "facts": {**_NONE, "has_barrel": True, "barrel_in_construct": 1,
                                                  "barrel_modelled": 1}},  # ChuA barrel
        "5ABC": {"p": 0.05, "why": [], "facts": {**_NONE, "has_barrel": True, "barrel_in_construct": 1,
                                                 "barrel_modelled": 1}}})
    v, ev, _, by = judge("7O7P", src)
    assert v == "uncertain" and by == "none" and "not enough" in ev[0]["text"]
    assert judge("9DIS", src)[:1] == ("membrane",) and judge("9DIS", src)[3] == "databases"
    assert judge("5ABC", src)[3] == "databases"  # resolved barrel: MemProtMD counts despite a low score


def test_swapped():
    assert vc.swapped("6IJ6", "6I6J") and vc.swapped("6SMX", "6XMS") and vc.swapped("7LIC", "7LCI")
    assert not vc.swapped("6IJ6", "6IJ6") and not vc.swapped("1ABC", "1BCA") and not vc.swapped("1ABC", "1ABD")
