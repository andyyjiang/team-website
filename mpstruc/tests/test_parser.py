import hashlib
from pathlib import Path

import pytest

from mpstruc import integrity
from mpstruc.flatten import structure_rows
from mpstruc.model import load, to_xml
from mpstruc.report import EXPECTED_SHA256
from mpstruc.resolution import parse_resolution
from mpstruc.text import mentioned_pdb_ids, to_text
from mpstruc.validate import compare, unexplained
from mpstruc.xmlparse import MpstrucParseError, parse, serialize

ROOT = Path(__file__).resolve().parent.parent
FIXTURE = ROOT / "fixtures" / "mpstruc_2026-10-02.xml"
XLSX = ROOT / "fixtures" / "mpstruc_XML.xlsx"


@pytest.fixture(scope="session")
def raw():
    return FIXTURE.read_bytes()


@pytest.fixture(scope="session")
def db(raw):
    return load(raw)


@pytest.fixture(scope="session")
def cmp(db):
    return compare(db, XLSX)


@pytest.fixture(scope="session")
def findings(db):
    return integrity.run(db, online=False)


def test_fixture_hash(raw):
    assert hashlib.sha256(raw).hexdigest() == EXPECTED_SHA256


def test_generic_round_trip_is_byte_identical(raw):
    assert serialize(parse(raw)).encode() == raw


def test_typed_round_trip_is_byte_identical(db, raw):
    assert to_xml(db).encode() == raw


def test_strict_xml_parser_rejects_the_source(raw):
    # Documents why we use our own reader: the export is not well-formed XML.
    import xml.etree.ElementTree as ET
    with pytest.raises(ET.ParseError):
        ET.fromstring(raw)


def test_totals(db):
    ents = list(db.entries())
    related = [r.pdb for e in ents for r in e.related]
    assert len(db.groups) == 3
    assert [len(g.subgroups) for g in db.groups] == [30, 17, 159]
    assert sum(e.kind == "parent" for e in ents) == 2070
    assert sum(e.kind == "daughter" for e in ents) == 2248
    assert len(related) == 6273
    assert len(ents) + len(related) == 10591
    assert len({e.pdbCode for e in ents} | set(related)) == 10552
    assert sum(1 + len(e.secondary) for e in ents) == 4331


def test_record_keys(db):
    keys = {e.key: e for e in db.entries()}
    assert keys["G01-S001-P0001"].pdbCode == "1PRH"
    assert keys["G01-S001-P0003-D001"].pdbCode == "6BL4"
    assert [r.pdb for r in keys["G01-S004-P0001-D001"].related] == ["1OJD", "1OJC", "1OJB", "1OJ9"]
    sg = {s.key: s.name for s in db.subgroups()}
    assert sg["G01-S005"] == sg["G03-S136"] == "Hydrolases"
    assert sg["G03-S140"] == "Oxidoreductases"


def test_daughters_point_at_their_parent(db):
    for e in db.entries():
        if e.kind == "parent":
            for d in e.daughters:
                assert d.masterProteinPdbCode == e.pdbCode, d.key


def test_totals_match_xlsx_source_sheet(cmp):
    assert all(t["match"] for t in cmp["totals"].values()), cmp["totals"]


def test_no_unexplained_xlsx_differences(cmp):
    assert unexplained(cmp) == []


def test_xlsx_defect_is_exactly_the_1prh_branch(cmp):
    s = cmp["structures"]
    assert s["xlsx_blank_rows"] == 16
    assert s["extra_in_xlsx"] == []
    assert list(cmp["structures_missing_by_parent"]) == ["G01-S001-P0001"]
    assert len(s["missing_from_xlsx"]) == 16


def test_references_match_xlsx_exactly(cmp):
    r = cmp["references"]
    assert (r["missing_from_xlsx"], r["extra_in_xlsx"], r["cells"]) == ([], [], [])


def test_cyclooxygenases_has_22_listings(db):
    rows = [r for r in structure_rows(db) if r["Subfamily"] == "Cyclooxygenases"]
    assert len(rows) == 22


@pytest.mark.parametrize("raw,value,kind", [
    ("3.5", 3.5, "number"),
    ("1.77 &Aring;", 1.77, "number"),
    ("1.8 &Aring; (P6<sub>3</sub> crystal form)", 1.8, "number"),
    ("3.4/3.7 &Aring;", 3.4, "multiple"),
    ("NMR Structure", None, "method"),
    ("NMR strucuture", None, "method"),
    ("solid-state NMR structure", None, "method"),
    ("by electron crystallography", None, "method"),
    ("-", None, "unparseable"),
    ("HEK293 cells", None, "unparseable"),
    (None, None, "missing"),
])
def test_resolution_grammar(raw, value, kind):
    r = parse_resolution(raw)
    assert (r.value, r.kind) == (value, kind)


def test_resolution_is_never_zero_for_missing(db):
    for e in db.entries():
        assert parse_resolution(e.resolution).value != 0


def test_text_conversion():
    assert to_text("2.8&nbsp;&Aring;: 3SQC <br /> P3<sub>2</sub>21") == "2.8 Å: 3SQC\nP3221"
    assert to_text("3.3 &Ariing;: 5KLS") == "3.3 Å: 5KLS"
    assert to_text("&alfa;-cyclodextrin") == "α-cyclodextrin"
    assert to_text("Yildiz &Ouml, K") == "Yildiz Ö, K"
    assert to_text("a<br /> <br />b") == "a\n\nb"


@pytest.mark.parametrize("raw,ids", [
    ("closed state, 10mM KCl, 2.33 &Aring;: 6W0I", ["6W0I"]),
    ("3.80&Aring;: 7vMS, 8Si3 &amp; 8Si4", ["7VMS", "8SI3", "8SI4"]),
    ("296A mutant, 3.12 Å: 8E6L", ["8E6L"]),
    ("Dark-adapted structure, 1.85 Å 6S63", ["6S63"]),
    ("+HQNO, 2.7 Å, 1KF6. +DNP-19, 3.6 Å, 1KFY", ["1KF6", "1KFY"]),
    ("3.30 Å7V2D", ["7V2D"]),
    ("2.36 Å; 8E5V", ["8E5V"]),
    ("2.15 Å: : 1PO0", ["1PO0"]),
    ("3.6 & Aring;:6QMB", ["6QMB"]),
    ("3.34 Aring;: 3AOC", ["3AOC"]),
    ("Former 3A2R superseded by 3VZT", []),
])
def test_mentioned_ids_need_resolution_context(raw, ids):
    assert mentioned_pdb_ids(raw) == ids


def test_text_repairs():
    assert to_text("x<br />\u200b\u200bsoaked") == "x\nsoaked"
    assert to_text("3.6 & Aring;:6QMB") == "3.6 Å:6QMB"


@pytest.mark.parametrize("old,new", [
    (b"<protein>\n<pdbCode>1PRH", b'<protein id="7">\n<pdbCode>1PRH'),
    (b"<description>In complex with bromoaspirin.</description>",
     b"<description>In complex with bromoaspirin.</description><description>x</description>"),
    (b"<pdbCode>1EQH</pdbCode>", b"<pdbCode/>"),
    (b"<pdbCode>1EQH</pdbCode>", b"<pdbCode> </pdbCode>"),
])
def test_parser_rejects_silent_corruptions(raw, old, new):
    assert raw.count(old) >= 1
    with pytest.raises(MpstrucParseError):
        load(raw.replace(old, new, 1))


def test_harmless_variants_are_normalised(raw):
    variant = b"\xef\xbb\xbf<?xml version=\"1.0\"?>\r\n" + raw.replace(b"\n", b"\r\n") + b"\r\n"
    db2 = load(variant)
    assert to_xml(db2).encode() == raw
    assert db2.normalizations


def test_parser_rejects_schema_violations(raw):
    bad = raw.replace(b"<masterProteinPdbCode>1PRH</masterProteinPdbCode>", b"", 1)
    with pytest.raises(MpstrucParseError):
        load(bad)


def _keys(findings, check):
    return {f["key"] for f in findings if f["check"] == check}


def test_integrity_cross_listing_kept(db, findings):
    occ = [f for f in findings if f["check"] == "cross_listing" and f["pdb"] == "3HYW"]
    assert occ and occ[0]["data"]["documented"]
    assert sum(1 for e in db.entries() if e.pdbCode == "3HYW") == 2


def test_integrity_known_issues(findings):
    assert {"G03-S038-P0017-D003", "G03-S038-P0017-D004"} <= _keys(findings, "malformed_markup")
    assert "G03-S039-P0004-D001" in _keys(findings, "unclosed_link")
    assert {"G03-S030-P0011", "G03-S037-P0018-D004"} <= _keys(findings, "resolution_unparseable")
    assert _keys(findings, "bad_pdb_id") == set()  # lower-case IDs are reported online as pdb_id_case
    assert "G03-S051-P0023" in _keys(findings, "described_not_related")  # cites 8J6D
    assert _keys(findings, "master_code_mismatch") == set()
    assert _keys(findings, "empty_subgroup") == {"G03-S050"}


def test_empty_header_subgroup_is_kept(db):
    sg = {s.key: s for s in db.subgroups()}
    assert sg["G03-S050"].name == "G PROTEIN-COUPLED RECEPTORS (GPCRs)" and sg["G03-S050"].parents == []
