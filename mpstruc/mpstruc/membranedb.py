"""Membrane-protein databases used as labels and signals.

- OPM (Orientations of Proteins in Membranes, opm.phar.umich.edu): every
  structure it positions in a membrane, with its type: transmembrane,
  monotopic/peripheral, or peptide. Read from OPM's own API, so it covers
  entries the PDB's copy of OPM no longer tags.
- TCDB (Transporter Classification Database, tcdb.org): UniProt accessions of
  classified transporters and channels, with their TC numbers.

PDBTM and MemProtMD are read through the PDB's search service (membrane.py).
"""
from __future__ import annotations

import json
import time
import urllib.request

from .rcsb import CACHE_DIR

OPM_API = "https://opm-back.cc.lehigh.edu/opm-backend"
TCDB_MAP = "https://www.tcdb.org/cgi-bin/projectv/public/acc2tcid.py"
OPM_TYPES = {1: "transmembrane", 2: "monotopic/peripheral", 3: "peptide"}


def _get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "mpstruc-curation-tool"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** (attempt + 1))
    raise RuntimeError("unreachable")


def _paged(endpoint: str, page_size: int = 1000) -> list[dict]:
    out, page = [], 1
    while True:
        d = json.loads(_get(f"{OPM_API}/{endpoint}?pageSize={page_size}&pageNum={page}"))
        d = d.get("table", d)
        out += d.get("objects") or []
        if page * page_size >= d.get("total_objects", 0):
            return out
        page += 1


def opm(max_age_s: float = 30 * 86400) -> dict[str, str]:
    """{PDB ID: 'transmembrane' | 'monotopic/peripheral' | 'peptide'} for every structure in OPM."""
    path = CACHE_DIR / "opm.json"
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        return json.loads(path.read_text())
    prim = _paged("primary_structures")
    by_id = {p["id"]: OPM_TYPES.get(p.get("type_id"), "other") for p in prim}
    out = {p["pdbid"].upper(): by_id[p["id"]] for p in prim if p.get("pdbid")}
    for s in _paged("secondary_representations"):
        if s.get("pdbid") and s.get("primary_structure_id") in by_id:
            out.setdefault(s["pdbid"].upper(), by_id[s["primary_structure_id"]])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out))
    return out


def tcdb(max_age_s: float = 30 * 86400) -> dict[str, list[str]]:
    """{UniProt accession: [TC numbers]} for every protein TCDB classifies."""
    path = CACHE_DIR / "tcdb.json"
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        return json.loads(path.read_text())
    out: dict[str, list[str]] = {}
    for line in _get(TCDB_MAP).decode("utf-8", "replace").splitlines():
        parts = line.strip().split("\t")
        if len(parts) == 2 and parts[0] and parts[1]:
            out.setdefault(parts[0], []).append(parts[1])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out))
    return out
