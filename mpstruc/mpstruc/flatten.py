"""Flatten a Database into the 4-sheet layout of the reference xlsx."""
from __future__ import annotations

from .model import Bibliography, Database, Entry
from .resolution import parse_resolution
from .text import to_text

STRUCTURE_COLUMNS = [
    "Subfamily", "PDB ID", "Topology", "Record type", "Protein name", "Description",
    "Species", "Taxonomic domain", "Expressed in species", "Resolution (Å)",
    "Resolution as recorded", "Parent PDB ID", "Related to PDB ID", "Record key",
    "Parent record key", "Related to record key", "Publication year", "DOI",
    "PubMed ID", "Related entry count", "Daughter record count", "PDB URL",
]
REFERENCE_COLUMNS = [
    "Record key", "PDB ID", "Subfamily", "Topology", "Reference type", "Reference number",
    "PubMed ID", "Authors", "Year", "Title", "Journal", "Volume", "Issue", "Pages", "DOI", "Notes",
]
SUBFAMILY_COLUMNS = [
    "Subfamily", "Topology", "Topology as recorded", "Parent records", "Daughter records",
    "Related PDB listings", "Total listings", "Distinct PDB IDs",
]


def _int(s: str | None) -> int | None:
    t = to_text(s)
    return int(t) if t and t.isdigit() else None


def _url(pdb: str) -> str:
    return f"https://www.rcsb.org/structure/{pdb}"


def structure_rows(db: Database) -> list[dict]:
    rows = []
    for g in db.groups:
        for sg in g.subgroups:
            for p in sg.parents:
                for e in [p, *p.daughters]:
                    res = parse_resolution(e.resolution)
                    b = e.bibliography
                    rows.append({
                        "Subfamily": to_text(sg.name), "PDB ID": e.pdbCode,
                        "Topology": g.topology,
                        "Record type": "Parent" if e.kind == "parent" else "Daughter",
                        "Protein name": to_text(e.name), "Description": to_text(e.description),
                        "Species": to_text(e.species), "Taxonomic domain": to_text(e.taxonomicDomain),
                        "Expressed in species": to_text(e.expressedInSpecies),
                        "Resolution (Å)": res.value, "Resolution as recorded": to_text(e.resolution),
                        "Parent PDB ID": p.pdbCode, "Related to PDB ID": None,
                        "Record key": e.key, "Parent record key": p.key, "Related to record key": None,
                        "Publication year": _int(b.year), "DOI": to_text(b.doi),
                        "PubMed ID": to_text(b.pubMedId), "Related entry count": len(e.related),
                        "Daughter record count": len(p.daughters) if e.kind == "parent" else None,
                        "PDB URL": _url(e.pdbCode),
                    })
                    for r in e.related:
                        rows.append({c: None for c in STRUCTURE_COLUMNS} | {
                            "Subfamily": to_text(sg.name), "PDB ID": r.pdb, "Topology": g.topology,
                            "Record type": "Related", "Parent PDB ID": p.pdbCode,
                            "Related to PDB ID": e.pdbCode, "Record key": r.key,
                            "Parent record key": p.key, "Related to record key": e.key,
                            "PDB URL": _url(r.pdb),
                        })
    return rows


def _ref_row(e: Entry, sg_name: str, topology: str, b: Bibliography, kind: str, n: int) -> dict:
    return {
        "Record key": e.key, "PDB ID": e.pdbCode, "Subfamily": sg_name, "Topology": topology,
        "Reference type": kind, "Reference number": n, "PubMed ID": to_text(b.pubMedId),
        "Authors": to_text(b.authors), "Year": _int(b.year), "Title": to_text(b.title),
        "Journal": to_text(b.journal), "Volume": to_text(b.volume), "Issue": to_text(b.issue),
        "Pages": to_text(b.pages), "DOI": to_text(b.doi), "Notes": to_text(b.notes),
    }


def reference_rows(db: Database) -> list[dict]:
    rows = []
    for g in db.groups:
        for sg in g.subgroups:
            for e in (x for p in sg.parents for x in (p, *p.daughters)):
                rows.append(_ref_row(e, to_text(sg.name), g.topology, e.bibliography, "Primary", 1))
                for i, b in enumerate(e.secondary, 2):  # numbered after the primary, as in the xlsx
                    rows.append(_ref_row(e, to_text(sg.name), g.topology, b, "Secondary", i))
    return rows


def subfamily_rows(db: Database) -> list[dict]:
    rows = []
    for g in db.groups:
        for sg in g.subgroups:
            ents = [x for p in sg.parents for x in (p, *p.daughters)]
            rel = [r.pdb for e in ents for r in e.related]
            ids = {e.pdbCode for e in ents} | set(rel)
            rows.append({
                "Subfamily": to_text(sg.name), "Subgroup key": sg.key, "Topology": g.topology,
                "Topology as recorded": g.name, "Parent records": len(sg.parents),
                "Daughter records": len(ents) - len(sg.parents), "Related PDB listings": len(rel),
                "Total listings": len(ents) + len(rel), "Distinct PDB IDs": len(ids),
            })
    return rows
