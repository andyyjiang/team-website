"""Compare our flattening of the XML with the reference xlsx, cell by cell.

Every difference is classified. Only "unexplained" differences indicate a
parser problem; the others are known defects or limitations of the xlsx.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

import openpyxl

from .flatten import reference_rows, structure_rows, subfamily_rows
from .model import Database
from .text import BROKEN_ARING_RE, ENTITY_REPAIRS, INVISIBLE_RE


def _sheet(wb, name: str) -> list[dict]:
    it = wb[name].iter_rows(values_only=True)
    header = next(it)
    return [dict(zip(header, r)) for r in it]


def _equal(ours, theirs) -> bool:
    if ours == theirs:
        return True
    if isinstance(theirs, (int, float)) and ours is not None:
        try:
            return float(ours) == float(theirs)
        except (TypeError, ValueError):
            return False
    return False


def classify(column: str, ours, theirs) -> str:
    if column == "Resolution (Å)" and theirs is None and ours is not None:
        return "resolution_parsed"  # xlsx left a parseable value such as '1.77 Å' blank
    if isinstance(theirs, str) and isinstance(ours, str):
        fixed = theirs
        for bad, good in ENTITY_REPAIRS.items():
            fixed = fixed.replace(f"&{bad};", good)
        fixed = BROKEN_ARING_RE.sub("Å", fixed)
        if fixed == ours:
            return "entity_repaired"  # xlsx kept a misspelled entity such as &Ariing; verbatim
        if INVISIBLE_RE.sub("", fixed) == ours:
            return "invisible_removed"  # xlsx kept zero-width characters
    return "unexplained"


def _diff_rows(ours: list[dict], theirs: list[dict], key) -> dict:
    o = {key(r): r for r in ours}
    t = {key(r): r for r in theirs if r.get("Record key")}
    blank_rows = sum(1 for r in theirs if not any(v is not None for v in r.values()))
    cells: list[dict] = []
    for k in o.keys() & t.keys():
        for col, tv in t[k].items():
            ov = o[k].get(col)
            if not _equal(ov, tv):
                cells.append({"key": k, "column": col, "ours": ov, "xlsx": tv,
                              "class": classify(col, ov, tv)})
    cells.sort(key=lambda c: (str(c["key"]), c["column"]))
    return {
        "ours": len(o), "xlsx_populated": len(t), "xlsx_blank_rows": blank_rows,
        "missing_from_xlsx": sorted(o.keys() - t.keys(), key=str),
        "extra_in_xlsx": sorted(t.keys() - o.keys(), key=str),
        "cells": cells,
        "cell_classes": dict(Counter(c["class"] for c in cells)),
    }


def compare(db: Database, xlsx_path: str | Path) -> dict:
    wb = openpyxl.load_workbook(xlsx_path, read_only=True)
    s_rows = structure_rows(db)
    structures = _diff_rows(s_rows, _sheet(wb, "Structures"), lambda r: r["Record key"])
    references = _diff_rows(
        reference_rows(db), _sheet(wb, "References"),
        lambda r: (r["Record key"], r["Reference type"], r["Reference number"]))

    # Subfamilies: matched by position, since names repeat (two "Hydrolases").
    ours_sf, theirs_sf = subfamily_rows(db), _sheet(wb, "Subfamilies")
    sf_cells = []
    for i, (o, t) in enumerate(zip(ours_sf, theirs_sf)):
        for col, tv in t.items():
            if not _equal(o.get(col), tv):
                sf_cells.append({"key": o["Subgroup key"], "column": col, "ours": o.get(col), "xlsx": tv})
    subfamilies = {"ours": len(ours_sf), "xlsx": len(theirs_sf), "cells": sf_cells}

    # Source sheet totals vs ours.
    source = {r["Field"]: r["Value"] for r in _sheet(wb, "Source") if r.get("Field")}
    ents = list(db.entries())
    related = [r.pdb for e in ents for r in e.related]
    totals = {
        "Parent records": sum(e.kind == "parent" for e in ents),
        "Daughter records": sum(e.kind == "daughter" for e in ents),
        "Related PDB listings": len(related),
        "All PDB listings": len(ents) + len(related),
        "Subfamilies": len(list(db.subgroups())),
        "Distinct parent PDB IDs": len({e.pdbCode for e in ents if e.kind == "parent"}),
        "Distinct PDB IDs in all listings": len({e.pdbCode for e in ents} | set(related)),
        "Bibliography records": sum(1 + len(e.secondary) for e in ents),
    }
    totals_cmp = {k: {"ours": v, "xlsx": source.get(k), "match": v == source.get(k)} for k, v in totals.items()}

    # What the rows missing from the Structures sheet are.
    by_key = {r["Record key"]: r for r in s_rows}
    missing_tree = defaultdict(list)
    for k in structures["missing_from_xlsx"]:
        missing_tree[by_key[k]["Parent record key"]].append(k)

    return {
        "xml_sha256": db.sha256,
        "xlsx_recorded_sha256": source.get("SHA-256"),
        "timestamp_match": db.meta.get("timeStamp") == source.get("Snapshot timestamp"),
        "totals": totals_cmp,
        "structures": structures,
        "structures_missing_by_parent": dict(missing_tree),
        "references": references,
        "subfamilies": subfamilies,
    }


def unexplained(report: dict) -> list[dict]:
    out = [c for s in ("structures", "references") for c in report[s]["cells"] if c["class"] == "unexplained"]
    out += report["subfamilies"]["cells"]
    return out
