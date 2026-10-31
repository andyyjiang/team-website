"""Deliverable 1 report: parser validation, xlsx diff and integrity findings.

    python -m mpstruc.report [--offline]

Writes reports/validation.md, reports/xlsx_diff.json, reports/integrity.json
and reports/integrity.md.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from . import integrity
from .model import load, to_xml
from .validate import compare, unexplained

ROOT = Path(__file__).resolve().parent.parent
FIXTURE = ROOT / "fixtures" / "mpstruc_2026-10-02.xml"
XLSX = ROOT / "fixtures" / "mpstruc_XML.xlsx"
EXPECTED_SHA256 = "4174f67bb5685afe0c0cf046f0c5bc5ceaa53162f8fd9568cf679f5187ea9f7a"


def source_xml() -> Path:
    """The mpstruc XML the website and automation use: $MPSTRUC_XML if set
    (the weekly job points it at a freshly downloaded copy), else the fixture."""
    import os
    p = os.environ.get("MPSTRUC_XML")
    return Path(p) if p and Path(p).exists() else FIXTURE

CHECK_TITLES = {
    "pdb_id_points_elsewhere": "PDB ID points to an unrelated structure",
    "pdb_title_mismatch": "PDB entry title shares no words with the protein name",
    "citation_ids_wrong": "Citation title matches the PDB's, but its PubMed ID/DOI do not",
    "citation_preprint": "mpstruc cites a preprint that has since been published",
    "pubmed_title_mismatch": "PubMed ID points to a different paper than the citation title",
    "pubmed_id_not_found": "PubMed IDs that do not exist",
    "pubmed_id_conflict": "DOI matches the PDB but the PubMed ID does not",
    "doi_conflict": "PubMed ID matches the PDB but the DOI does not",
    "resolution_vs_method": "Resolution conflicts with the PDB's experimental method",
    "broken_entity": "Malformed entity spellings",
    "invisible_characters": "Zero-width characters",
    "citation_differs_from_pdb": "Primary citation differs from the PDB entry's citation",
    "resolution_differs_from_pdb": "Resolution differs from the PDB entry",
    "pdb_obsolete": "Obsolete or superseded PDB IDs",
    "pdb_unknown": "PDB IDs the PDB does not recognise",
    "pdb_id_case": "PDB IDs in the wrong case",
    "pdb_unreleased": "PDB IDs not yet released",
    "bad_pdb_id": "Malformed PDB IDs",
    "master_code_mismatch": "Daughter's masterProteinPdbCode differs from its parent",
    "self_related": "Entries listing themselves as related",
    "repeated_in_subgroup": "PDB IDs listed more than once in the same subgroup",
    "cross_listing": "PDB IDs listed in more than one subgroup (usually intentional)",
    "described_not_related": "Description cites a structure missing from the related list",
    "related_not_described": "Related entries not cited in the description",
    "duplicate_related": "Repeated related entries",
    "resolution_unparseable": "Resolution text that is neither a number nor a method",
    "resolution_missing": "Missing resolution",
    "resolution_multiple": "Resolution with several values",
    "unknown_domain": "Taxonomic domain outside Eukaryota/Bacteria/Archaea/Viruses",
    "missing_field": "Empty name, species or taxonomic domain",
    "citation_incomplete": "Citation without a title",
    "citation_unlinked": "Citation with neither PubMed ID nor DOI",
    "bad_pubmed_id": "Malformed PubMed IDs",
    "bad_doi": "Malformed DOIs",
    "bad_year": "Implausible publication years",
    "unknown_entity": "Misspelled HTML entities",
    "unterminated_entity": "HTML entities missing ';'",
    "bare_ampersand": "Unescaped '&' (makes the file invalid XML)",
    "malformed_markup": "Unexpected or malformed tags",
    "unclosed_link": "Unclosed <a> links",
    "empty_tag": "Empty formatting tags",
    "empty_subgroup": "Subgroups with no entries",
    "duplicate_subgroup_name": "Subgroup names used more than once",
}
SEVERITY_ORDER = {"error": 0, "warning": 1, "info": 2}


def _fmt(v) -> str:
    s = "" if v is None else str(v)
    s = s.replace("\n", " ⏎ ").replace("|", "\\|")
    return s if len(s) <= 140 else s[:137] + "…"


def integrity_md(findings: list[dict], db) -> str:
    by = defaultdict(list)
    for f in findings:
        by[f["check"]].append(f)
    order = sorted(by, key=lambda c: (min(SEVERITY_ORDER[f["severity"]] for f in by[c]), c))
    lines = [
        "# mpstruc data-integrity report", "",
        f"Snapshot `{db.meta.get('timeStamp')}`, SHA-256 `{db.sha256}`.", "",
        "Findings are for curator review. Nothing in the data has been changed.", "",
        "| Check | Severity | Count |", "|---|---|---|",
    ]
    for c in order:
        sev = Counter(f["severity"] for f in by[c])
        lines.append(f"| [{CHECK_TITLES.get(c, c)}](#{c}) | "
                     + ", ".join(f"{k} {v}" for k, v in sorted(sev.items(), key=lambda kv: SEVERITY_ORDER[kv[0]]))
                     + f" | {len(by[c])} |")
    for c in order:
        lines += ["", f"<a id=\"{c}\"></a>", f"## {CHECK_TITLES.get(c, c)}", "",
                  "| Severity | Record key | PDB | Finding |", "|---|---|---|---|"]
        for f in sorted(by[c], key=lambda f: (SEVERITY_ORDER[f["severity"]], f["key"] or "")):
            lines.append(f"| {f['severity']} | `{f['key']}` | {f['pdb'] or ''} | {_fmt(f['message'])} |")
    return "\n".join(lines) + "\n"


def validation_md(cmp: dict, roundtrip_ok: bool, findings: list[dict], online: bool) -> str:
    s, r = cmp["structures"], cmp["references"]
    unexp = unexplained(cmp)
    sev = Counter(f["severity"] for f in findings)
    L = [
        "# Deliverable 1: parser validation", "",
        "## Summary", "",
        f"- **Round trip:** parsing the fixture and writing it back is **{'byte-identical' if roundtrip_ok else 'NOT identical'}**.",
        f"- **Totals:** {'all' if all(t['match'] for t in cmp['totals'].values()) else 'NOT all'} "
        "match the xlsx `Source` sheet (table below).",
        f"- **Cell-by-cell diff against the xlsx:** {len(unexp)} unexplained difference(s). "
        f"{sum(len(x['cells']) for x in (s, r)) - len(unexp) + len(cmp['subfamilies']['cells'])} explained ones are listed below.",
        f"- **xlsx defect:** the `Structures` sheet has {s['xlsx_blank_rows']} blank rows; the "
        f"{len(s['missing_from_xlsx'])} records missing from it are the whole "
        f"{', '.join(cmp['structures_missing_by_parent'])} branch (1PRH, ram COX-1). They are present in the XML.",
        f"- **Integrity findings:** {sev.get('error', 0)} errors, {sev.get('warning', 0)} warnings, "
        f"{sev.get('info', 0)} info (see `integrity.md`)"
        + ("" if online else "; PDB cross-checks were skipped (offline run)") + ".",
        "",
        "## Snapshot identity", "",
        f"- Fixture SHA-256: `{cmp['xml_sha256']}`",
        f"- SHA-256 recorded in the xlsx `Source` sheet: `{cmp['xlsx_recorded_sha256']}`. These differ. The "
        f"timestamps {'match' if cmp['timestamp_match'] else 'DO NOT match'} and every count below matches, so the "
        "xlsx hash was most likely taken over different bytes (for example before re-saving). The fixture hash above is "
        "the one the parser checks.",
        "",
        "## Totals", "",
        "| Quantity | Parsed from XML | xlsx `Source` sheet | Match |", "|---|---|---|---|",
    ]
    for k, v in cmp["totals"].items():
        L.append(f"| {k} | {v['ours']:,} | {v['xlsx']:,} | {'✓' if v['match'] else '✗'} |")
    L += [
        "", "## Structures sheet", "",
        f"- Parsed listings: {s['ours']:,}; populated xlsx rows: {s['xlsx_populated']:,}; blank xlsx rows: {s['xlsx_blank_rows']}",
        f"- Rows in the xlsx that we do not produce: {len(s['extra_in_xlsx'])}",
        f"- Records missing from the xlsx: {len(s['missing_from_xlsx'])}",
        "",
        "Missing records, grouped by parent:", "",
    ]
    for p, ks in cmp["structures_missing_by_parent"].items():
        L.append(f"- `{p}`: {', '.join('`' + k + '`' for k in ks)}")
    L += ["", "### Explained cell differences", "",
          "| Class | Meaning | Count |", "|---|---|---|",
          f"| resolution_parsed | The xlsx left a parseable resolution blank (e.g. `1.77 Å`, `1.8 Å (P63 crystal form)`, `3.4/3.7 Å`) | {s['cell_classes'].get('resolution_parsed', 0)} |",
          f"| entity_repaired | The xlsx kept a misspelled entity (`&Ariing;`, `&alfa;`, `& Aring;`, `Aring;`) as literal text; we display Å / α | {s['cell_classes'].get('entity_repaired', 0)} |",
          f"| invisible_removed | The xlsx kept zero-width characters | {s['cell_classes'].get('invisible_removed', 0)} |",
          "", "| Record key | Column | Parsed | xlsx | Class |", "|---|---|---|---|---|"]
    for c in s["cells"]:
        L.append(f"| `{c['key']}` | {c['column']} | {_fmt(c['ours'])} | {_fmt(c['xlsx'])} | {c['class']} |")
    L += ["", "## References sheet", "",
          f"- Parsed bibliography records: {r['ours']:,}; xlsx rows: {r['xlsx_populated']:,}",
          f"- Missing / extra / differing: {len(r['missing_from_xlsx'])} / {len(r['extra_in_xlsx'])} / {len(r['cells'])}",
          "- The 9 reference rows for the 1PRH branch are present in this sheet even though the branch's "
          "structure rows are missing from `Structures`.",
          "", "## Subfamilies sheet", "",
          f"- Parsed subgroups: {cmp['subfamilies']['ours']}; xlsx rows: {cmp['subfamilies']['xlsx']}; "
          f"differing cells: {len(cmp['subfamilies']['cells'])}",
          "- Subgroups are numbered within their group (`G03-S136`), matching the xlsx record keys. Two subgroups are "
          "named *Hydrolases* (G01-S005, G03-S136), so subgroups are identified by key, never by name.",
          "", "## Known limitations of the xlsx (kept in our data, lost in the xlsx)", "",
          "- Formatting is stripped: `P4<sub>3</sub>2<sub>1</sub>2` becomes `P43212`, and `Mn<sup>2+</sup>` becomes `Mn2+`. "
          "We keep the raw markup and render it.",
          "- Related entries carry no metadata in the XML. We will fill it in from the PDB at display time, labelled \"from PDB\".",
          ""]
    return "\n".join(L)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true", help="skip RCSB cross-checks")
    ap.add_argument("--out", default=str(ROOT / "reports"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    raw = FIXTURE.read_bytes()
    db = load(raw)
    if db.sha256 != EXPECTED_SHA256:
        raise SystemExit(f"fixture hash {db.sha256} != expected {EXPECTED_SHA256}")
    roundtrip_ok = to_xml(db).encode("utf-8") == raw
    cmp = compare(db, XLSX)
    findings = integrity.run(db, online=not args.offline)

    (out / "xlsx_diff.json").write_text(json.dumps(cmp, indent=1, default=str, ensure_ascii=False))
    (out / "integrity.json").write_text(json.dumps(findings, indent=1, ensure_ascii=False))
    (out / "integrity.md").write_text(integrity_md(findings, db))
    (out / "validation.md").write_text(validation_md(cmp, roundtrip_ok, findings, not args.offline))
    print(f"round trip {'OK' if roundtrip_ok else 'FAILED'}; unexplained xlsx diffs: {len(unexplained(cmp))}; "
          f"findings: {dict(Counter(f['severity'] for f in findings))}; reports in {out}")


if __name__ == "__main__":
    main()
