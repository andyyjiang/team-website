"""How do mpstruc curators place structures? Measured from the snapshot.

    python -m mpstruc.placement      # writes reports/placement_analysis.{md,json}

For every daughter we compare it with its parent, and for every related
entry with its owner (the parent or daughter it is listed under): do they
share a membrane-protein UniProt accession, the same organism, the same
primary citation? The answers set the placement rules the automation uses.

"Membrane entities" are the polymer entities the PDB annotates as membrane
proteins (OPM, PDBTM or MemProtMD). RCSB's copy of the mpstruc
classification is deliberately not used, so that the rules can be tested
on entries the automation has not seen.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from .model import Database, load
from .report import FIXTURE, ROOT

MEMBRANE_SOURCES = {"OPM", "PDBTM", "MemProtMD"}


def entity_summary(entry: dict | None) -> dict | None:
    """UniProt IDs, organisms, Pfam IDs and citation for an entry's membrane entities."""
    if not entry:
        return None
    ents = entry.get("polymer_entities") or []
    proteins = [e for e in ents if ((e.get("entity_poly") or {}).get("rcsb_entity_polymer_type") == "Protein")]
    membrane = [e for e in proteins
                if any(a.get("type") in MEMBRANE_SOURCES for a in e.get("rcsb_polymer_entity_annotation") or [])]
    chosen = membrane or proteins

    def ids(es, f):
        out = set()
        for e in es:
            out |= set(f(e) or [])
        return out

    cit = entry.get("rcsb_primary_citation") or {}
    pm = cit.get("pdbx_database_id_PubMed")
    return {
        "uniprot": ids(chosen, lambda e: (e.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids")),
        "all_uniprot": ids(proteins, lambda e: (e.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids")),
        "taxa": ids(chosen, lambda e: [o.get("ncbi_taxonomy_id") for o in e.get("rcsb_entity_source_organism") or []
                                       if o.get("ncbi_taxonomy_id")]),
        "pfam": ids(chosen, lambda e: [a["annotation_id"] for a in e.get("rcsb_polymer_entity_annotation") or []
                                       if a.get("type") == "Pfam"]),
        "membrane_annotated": bool(membrane),
        "sources": ids(proteins, lambda e: [a["type"] for a in e.get("rcsb_polymer_entity_annotation") or []
                                            if a.get("type") in MEMBRANE_SOURCES]),
        "citation": {x for x in (str(cit.get("pdbx_database_id_DOI") or "").lower(),
                                 str(pm) if pm and int(pm) > 0 else "") if x},
        "released": (entry.get("rcsb_accession_info") or {}).get("initial_release_date", "")[:10],
        "n_protein_entities": len(proteins),
    }


def relation(a: dict | None, b: dict | None) -> dict | None:
    """How structure a relates to structure b."""
    if not a or not b:
        return None
    return {
        "same_uniprot": bool(a["uniprot"] & b["uniprot"]),
        "any_shared_uniprot": bool(a["all_uniprot"] & b["all_uniprot"]),
        "same_taxon": bool(a["taxa"] & b["taxa"]),
        "same_citation": bool(a["citation"] & b["citation"]),
        "same_pfam": bool(a["pfam"] & b["pfam"]),
        "same_uniprot_set": bool(a["uniprot"]) and a["uniprot"] == b["uniprot"],
    }


def analyse(db: Database, entities: dict[str, dict | None]) -> dict:
    S = {k.upper(): entity_summary(v) for k, v in entities.items()}
    rows = {"daughter_vs_parent": [], "related_vs_owner": []}
    for sg in db.subgroups():
        for p in sg.parents:
            sp = S.get(p.pdbCode.upper())
            for d in p.daughters:
                rows["daughter_vs_parent"].append(relation(S.get(d.pdbCode.upper()), sp))
            for owner in (p, *p.daughters):
                so = S.get(owner.pdbCode.upper())
                for r in owner.related:
                    rows["related_vs_owner"].append(relation(S.get(r.pdb.upper()), so))

    def tally(rs):
        rs = [r for r in rs if r]
        n = len(rs)
        c = Counter()
        for r in rs:
            for k, v in r.items():
                c[k] += v
            c["uniprot+citation"] += r["same_uniprot"] and r["same_citation"]
            c["uniprot, other citation"] += r["same_uniprot"] and not r["same_citation"]
            c["other uniprot, same taxon"] += (not r["same_uniprot"]) and r["same_taxon"]
            c["other uniprot, other taxon"] += (not r["same_uniprot"]) and not r["same_taxon"]
        return {"n": n, **{k: {"count": v, "pct": round(100 * v / n, 1) if n else 0} for k, v in c.items()}}

    # Parents: is a new parent ever the same protein as an earlier parent?
    parents = [(p, S.get(p.pdbCode.upper())) for sg in db.subgroups() for p in sg.parents]
    by_uniprot: dict[str, list] = {}
    for p, s in parents:
        for u in (s or {}).get("uniprot", ()):
            by_uniprot.setdefault(u, []).append(p.key)
    shared_parent_uniprot = {u: ks for u, ks in by_uniprot.items() if len(ks) > 1}

    coverage = Counter()
    for k, s in S.items():
        coverage["entries"] += 1
        coverage["found in PDB"] += s is not None
        coverage["membrane-annotated (OPM/PDBTM/MemProtMD)"] += bool(s and s["membrane_annotated"])
        coverage["has UniProt for membrane entities"] += bool(s and s["uniprot"])

    return {
        "daughter_vs_parent": tally(rows["daughter_vs_parent"]),
        "related_vs_owner": tally(rows["related_vs_owner"]),
        "parents_sharing_uniprot": {"uniprot_ids": len(shared_parent_uniprot),
                                    "parents": sum(len(v) for v in shared_parent_uniprot.values()),
                                    "examples": dict(list(shared_parent_uniprot.items())[:15])},
        "coverage": dict(coverage),
    }


def to_md(a: dict) -> str:
    L = ["# How mpstruc curators place structures", "",
         "Measured on the 2 Oct 2026 snapshot using PDB data (UniProt, organism, citation). "
         "'Same protein' means the two structures share a UniProt accession on their membrane entities.", ""]
    for key, title in (("daughter_vs_parent", "Daughters compared with their parent"),
                       ("related_vs_owner", "Related entries compared with the entry they are listed under")):
        t = a[key]
        L += [f"## {title} (n = {t['n']:,})", "", "| Relationship | Count | % |", "|---|---|---|"]
        for k in ("same_uniprot", "same_taxon", "same_citation", "same_pfam", "uniprot+citation",
                  "uniprot, other citation", "other uniprot, same taxon", "other uniprot, other taxon"):
            if k in t:
                L.append(f"| {k.replace('_', ' ')} | {t[k]['count']:,} | {t[k]['pct']} |")
        L.append("")
    p = a["parents_sharing_uniprot"]
    L += ["## Parents that are the same protein as another parent", "",
          f"{p['uniprot_ids']} UniProt accessions are shared by {p['parents']} parent records "
          "(e.g. cross-listings, or the same protein given a new parent for a new state or complex).", "",
          "## Coverage", ""] + [f"- {k}: {v:,}" for k, v in a["coverage"].items()]
    return "\n".join(L) + "\n"


def main() -> None:
    from . import rcsb
    db = load(FIXTURE)
    ids = sorted({e.pdbCode.upper() for e in db.entries()} | {r.pdb.upper() for e in db.entries() for r in e.related})
    a = analyse(db, rcsb.entry_entities(ids))
    out = ROOT / "reports"
    (out / "placement_analysis.json").write_text(json.dumps(a, indent=1, default=sorted))
    (out / "placement_analysis.md").write_text(to_md(a))
    print(to_md(a))


if __name__ == "__main__":
    main()
