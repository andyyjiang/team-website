import json
import shutil
import subprocess
from pathlib import Path

import pytest

from mpstruc import vercel_import as vi

SAMPLE = Path(__file__).resolve().parent.parent / "vercel" / "sample_page.txt"
READER = SAMPLE.parent / "read_page.js"

CSV = """﻿Pdb Code,Name,Experimental Method,Species,Taxonomic Domain,Expressed In Species,Resolution,Description,Family,Subfamily,Bibliography,PDB Link,mpstruc Link,PubMed Link,DOI Link,Status,Memo
9tgs,Cryo-EM structure of TamA,Cryo-EM,,,,7.75,,Outer-Membrane,"Type II,,Omp85",Paper.. Cell Rep (2026). 45:1. doi:10.1/x. PMID:42507553,https://www.rcsb.org/structure/9TGS,,https://www.ncbi.nlm.nih.gov/pubmed/?term=42507553,,,
1ABC,Some channel,X-ray,,,,2.1,,Multi-Pass,,No pmid here,,,https://www.ncbi.nlm.nih.gov/pubmed/?term=111,,Already in,
2ABC,"A soluble, thing",X-ray,,,,2.1,,non-membrane,,Ref. PMID:222,,,,,No,"soluble domain only"
3ABC,Ready one,X-ray,,,,2.1,,,,Ref. PMID:333,,,,,Ready for yes,
4ABC,Pub one,X-ray,,,,2.1,,Single-Pass,,Ref. PMID:444,,,,,Pubmed ready,
"""


@pytest.fixture
def export(tmp_path):
    p = tmp_path / "mpstruc_form_export_2026-11-13.csv"
    p.write_text(CSV, encoding="utf-8")
    return p


@pytest.fixture
def offline_dates(monkeypatch, tmp_path):
    dates = {"9TGS": "2026-11-11", "1ABC": "2001-01-01", "2ABC": "2020-02-02", "3ABC": "2026-11-04",
             "4ABC": "2026-11-04", "9TYV": "2026-10-07"}

    def fake(ids, batch=300, offline=False):
        return {p.upper(): ({"released": dates[p.upper()], "title": "RCSB " + p.upper()}
                            if p.upper() in dates else None) for p in ids}
    monkeypatch.setattr(vi, "release_info", fake)


def test_status_codes():
    assert vi.status_code("Already in") == "already_in"
    assert vi.status_code(" Pubmed ready ") == "pubmed_ready"
    assert vi.status_code("Ready for yes") == "ready_for_yes"
    assert vi.status_code("") == "pending"
    assert vi.status_code(None) == "pending"
    assert vi.status_code("Maybe") == "maybe"
    assert vi.status_code("Something New") == "something_new"


def test_read_export(export):
    rows = vi.read_export(export)
    assert set(rows) == {"9TGS", "1ABC", "2ABC", "3ABC", "4ABC"}
    assert rows["9TGS"] == {"status": "pending", "family": "Outer-Membrane", "name": "Cryo-EM structure of TamA",
                            "group": "42507553", "memo": None}
    assert rows["1ABC"]["group"] == "111"          # falls back to the PubMed link
    assert rows["2ABC"]["memo"] == "soluble domain only"
    assert rows["2ABC"]["name"] == "A soluble, thing"
    assert rows["3ABC"]["family"] is None


def test_import_export_fills_dates_and_keeps_page_fields(export, offline_dates):
    snap = vi.load_snapshot(Path("/nonexistent"))
    snap["entries"]["9TGS"] = {"status": "yes", "family": None, "name": None, "group": None, "released": None,
                               "tm": "N/A", "memo": None, "seenAt": "2026-11-12"}
    res = vi.import_export(export, snap, offline=True)
    assert res["rows"] == 5 and res["noReleaseDate"] == []
    e = snap["entries"]["9TGS"]
    assert e["status"] == "pending"                 # the export replaces the fields it has
    assert e["tm"] == "N/A"                         # ...and keeps the ones it lacks
    assert e["seenAt"] == "2026-11-12"              # first sighting kept
    assert e["released"] == "2026-11-11"
    assert snap["entries"]["1ABC"]["seenAt"] == "2026-11-13"   # export date from the file name
    assert snap["source"] == "export:mpstruc_form_export_2026-11-13.csv"
    assert set(snap["entries"]["1ABC"]) == {"status", "family", "name", "group", "released", "tm", "memo", "seenAt"}


def test_parse_page_text_sample():
    groups = vi.parse_page_text(SAMPLE.read_text(encoding="utf-8"))
    assert len(groups) == 50
    g = groups[0]
    assert g["group"] == "42507553"
    assert g["journal"] == "(2026) Cell Rep 45: 117733"
    assert g["paperTitle"] == "Mechanistic insights into TAM-mediated OMP assembly in Gram-negative bacteria."
    assert g["doi"] == "https://doi.org/10.1016/j.celrep.2026.117733"
    assert g["title"] == "Cryo-EM structure of TamA"
    assert g["pdbIds"] == ["9TGS", "9TGR", "9TGQ", "9TGP"]
    assert g["releaseDates"] == ["2026-11-11"]
    assert g["tm"] == "N/A" and g["resolution"] == "7.75 Å"
    assert g["uniprot"] == ["P0ADE4", "P39321"]
    assert g["classifications"] == ["Outer-Membrane"]
    assert g["predictions"][2] == {"model": "Text model on abstract",
                                   "label": "Omp85-TpsB Outer Membrane Transporter Superfamily", "pct": 41.13}
    assert g["verdict"] == "Yes"
    ac9 = next(x for x in groups if x["group"] == "42752937")
    assert ac9["releaseDates"] == ["2026-11-11", "2026-11-04"]
    assert ac9["tm"] == "0 (some N/A)"
    assert ac9["classifications"] == ["Multi-Pass", "non-membrane"]
    assert ac9["verdict"] == "(Pending)"
    assert {x["verdict"] for x in groups} == {"Yes", "No", "(Pending)"}
    assert all(x["pdbIds"] and x["releaseDates"] for x in groups)
    # the "Load 50 more groups" line after the last card is not taken for data
    assert groups[-1]["pdbIds"] == ["9TYV"] and groups[-1]["verdict"] == "Yes"


def test_merge_page(offline_dates):
    snap = vi.load_snapshot(Path("/nonexistent"))
    snap["entries"]["9TGS"] = {"status": "ready_for_yes", "family": "Outer-Membrane", "name": "x", "group": None,
                               "released": "2026-11-11", "tm": None, "memo": "m", "seenAt": "2026-11-13"}
    snap["entries"]["9TGR"] = {"status": "pending", "family": None, "name": None, "group": None,
                               "released": None, "tm": None, "memo": None, "seenAt": "2026-11-13"}
    snap["entries"]["30BB"] = {"status": "already_in", "family": None, "name": None, "group": None,
                               "released": None, "tm": None, "memo": None, "seenAt": "2026-11-13"}
    snap["entries"]["11LW"] = {"status": "yes", "family": None, "name": None, "group": None,
                               "released": None, "tm": None, "memo": None, "seenAt": "2026-11-13"}
    groups = vi.parse_page_text(SAMPLE.read_text(encoding="utf-8"))
    res = vi.merge_page(groups, snap, seen_day="2026-11-12", offline=True)
    e = snap["entries"]
    assert e["9TGS"]["status"] == "ready_for_yes"   # page "Yes" does not undo a more specific export status
    assert e["9TGR"]["status"] == "yes"
    assert e["30BB"]["status"] == "already_in"      # page "(Pending)" never overwrites a decided status
    assert e["11LW"]["status"] == "no"              # the page's newer No wins over Yes
    assert e["9TGS"]["memo"] == "m" and e["9TGS"]["group"] == "42507553" and e["9TGS"]["tm"] == "N/A"
    assert e["9TGS"]["seenAt"] == "2026-11-12"
    assert e["9TYV"]["released"] == "2026-10-07" and e["9TYV"]["name"] == "RCSB 9TYV"
    # a group with one release date fills dates RCSB does not know
    assert e["9TGQ"]["released"] == "2026-11-11"
    assert res["groups"] == 50 and res["added"] == res["pdbs"] - 4
    assert snap["source"] == "page" and snap["pageReadAt"]


def test_load_page_file_jsonl(tmp_path):
    groups = vi.parse_page_text(SAMPLE.read_text(encoding="utf-8"))[:2]
    p = tmp_path / "g.jsonl"
    p.write_text("".join(json.dumps(g) + "\n" for g in groups)
                 + json.dumps({"_meta": {"readAt": "2026-11-13T00:00:00Z", "total": 7098}}) + "\n")
    got, meta = vi.load_page_file(p)
    assert got == groups and meta["total"] == 7098
    got, meta = vi.load_page_file(SAMPLE)
    assert len(got) == 50 and meta == {}


@pytest.mark.skipif(not shutil.which("node"), reason="node not installed")
def test_js_reader_parses_like_python():
    out = subprocess.run(["node", str(READER), "--parse-text", str(SAMPLE)], capture_output=True, text=True,
                         check=True).stdout
    lines = [json.loads(ln) for ln in out.splitlines() if ln.strip()]
    assert lines[-1]["_meta"]["groups"] == 50
    assert lines[:-1] == vi.parse_page_text(SAMPLE.read_text(encoding="utf-8"))


def test_merge_page_already_in_verdict(offline_dates):
    snap = vi.load_snapshot(Path("/nonexistent"))
    g = {"group": "1", "pdbIds": ["5ABC"], "releaseDates": ["2026-11-04"], "verdict": "Already in",
         "classifications": ["Multi-Pass"], "tm": "12", "title": "T"}
    vi.merge_page([g], snap, seen_day="2026-11-13", offline=True)
    e = snap["entries"]["5ABC"]
    assert e["status"] == "already_in" and e["family"] == "Multi-Pass" and e["released"] == "2026-11-04"
    assert e["name"] == "T" and e["seenAt"] == "2026-11-13"


def test_empty_page_read_changes_nothing(tmp_path):
    path = tmp_path / "snap.json"
    before = {"source": "export:x.csv", "readAt": "2026-11-13T00:00:00+00:00", "export": "x.csv", "pageReadAt": None,
              "entries": {"1ABC": {"status": "no", "family": None, "name": None, "group": None, "released": None,
                                   "tm": None, "memo": None, "seenAt": "2026-11-13"}}}
    path.write_text(json.dumps(before))
    page = tmp_path / "page.jsonl"
    page.write_text(json.dumps({"_meta": {"groups": 0, "clicks": 2, "stopReason": "end of list"}}) + "\n")
    with pytest.raises(SystemExit) as e:
        vi.main(["page", str(page), "--offline", "--snapshot", str(path)])
    assert e.value.code not in (0, None)
    assert json.loads(path.read_text()) == before
    with pytest.raises(ValueError):
        vi.merge_page([{"group": "1", "pdbIds": []}], vi.load_snapshot(path), offline=True)


def test_memos_are_not_published(tmp_path):
    path, memos = tmp_path / "snap.json", tmp_path / "memos.json"
    snap = {"source": "export:x.csv", "readAt": None, "export": "x.csv", "pageReadAt": None,
            "entries": {"20XK": {"status": "no", "family": None, "name": "N", "group": None, "released": None,
                                 "tm": None, "memo": "not a membrane protein", "seenAt": "2026-11-13"}}}
    vi.save_snapshot(snap, path, memos=memos)
    assert "not a membrane protein" not in path.read_text()
    assert json.loads(path.read_text())["entries"]["20XK"]["memo"] is None
    assert vi.load_snapshot(path, memos=memos)["entries"]["20XK"]["memo"] == "not a membrane protein"


def test_parse_group_skips_lines_above_the_label():
    g = vi.parse_group("#1\nPubMed Group: 41\nJournal:\nNature\nPaper\nTitle: X\nPDB IDs:\n9AB1\n"
                       "Release Date(s): 2026-11-11\nMembrane Protein Verdict: No")
    assert g["group"] == "41" and g["pdbIds"] == ["9AB1"] and g["verdict"] == "No"
