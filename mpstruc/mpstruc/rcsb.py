"""Minimal cached RCSB PDB client (holdings only, for now)."""
from __future__ import annotations

import http.client
import json
import os
import time
import urllib.error
import urllib.request
from pathlib import Path

CACHE_DIR = Path(os.environ.get("MPSTRUC_CACHE", Path(__file__).resolve().parent.parent / "cache"))
BASE = "https://data.rcsb.org/rest/v1/holdings"


def _get_json(url: str, cache_name: str, max_age_s: float) -> object | None:
    path = CACHE_DIR / cache_name
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        return json.loads(path.read_text())
    req = urllib.request.Request(url, headers={"User-Agent": "mpstruc-curation-tool"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.load(r)
            break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                data = None
                break
            if attempt == 3:
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError, http.client.HTTPException):
            if attempt == 3:
                raise
        time.sleep(2 ** (attempt + 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return data


def current_entry_ids(max_age_s: float = 86400) -> set[str]:
    return set(_get_json(f"{BASE}/current/entry_ids", "current_entry_ids.json", max_age_s) or [])


def unreleased_entry_ids(max_age_s: float = 86400) -> set[str]:
    return set(_get_json(f"{BASE}/unreleased/entry_ids", "unreleased_entry_ids.json", max_age_s) or [])


def removed_entry(pdb_id: str, max_age_s: float = 30 * 86400) -> dict | None:
    """Holdings record for an obsoleted entry, or None if RCSB has none."""
    data = _get_json(f"{BASE}/removed/{pdb_id}", f"removed/{pdb_id}.json", max_age_s)
    return (data or {}).get("rcsb_repository_holdings_removed")


GRAPHQL = "https://data.rcsb.org/graphql"
_SUMMARY_FIELDS = """rcsb_id
  struct { title }
  exptl { method }
  rcsb_entry_info { resolution_combined }
  rcsb_accession_info { initial_release_date }
  rcsb_primary_citation { pdbx_database_id_DOI pdbx_database_id_PubMed title year journal_abbrev }"""


def _post_json(url: str, payload: dict) -> dict:
    body = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=body, headers={
        "Content-Type": "application/json", "User-Agent": "mpstruc-curation-tool"})
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except (urllib.error.URLError, TimeoutError, ConnectionError, http.client.HTTPException):
            if attempt == 5:
                raise
            time.sleep(2 ** (attempt + 1))
    raise RuntimeError("unreachable")


def entry_summaries(pdb_ids: list[str], batch: int = 200, max_age_s: float = 7 * 86400) -> dict[str, dict]:
    """Title, method, resolution, release date and primary citation per entry.

    Cached as one JSON file; only IDs missing from the cache are fetched.
    Entries the PDB does not know (obsolete, typos) are cached as None.
    """
    path = CACHE_DIR / "entry_summaries.json"
    cache: dict[str, dict | None] = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    todo = sorted({p.upper() for p in pdb_ids} - cache.keys())
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        q = "query($ids: [String!]!) { entries(entry_ids: $ids) { %s } }" % _SUMMARY_FIELDS
        data = _post_json(GRAPHQL, {"query": q, "variables": {"ids": chunk}})
        found = {e["rcsb_id"]: e for e in (data.get("data") or {}).get("entries") or [] if e}
        for pid in chunk:
            cache[pid] = found.get(pid)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cache))
    return {p: cache.get(p.upper()) for p in pdb_ids}


_ENTITY_FIELDS = """rcsb_id
  rcsb_accession_info { initial_release_date deposit_date }
  rcsb_primary_citation { pdbx_database_id_DOI pdbx_database_id_PubMed }
  polymer_entities {
    rcsb_id
    rcsb_polymer_entity { pdbx_description }
    rcsb_polymer_entity_container_identifiers { uniprot_ids auth_asym_ids }
    rcsb_entity_source_organism { ncbi_taxonomy_id ncbi_scientific_name }
    entity_poly { rcsb_entity_polymer_type }
    rcsb_polymer_entity_annotation { type annotation_id name }
  }"""


def entry_entities(pdb_ids: list[str], batch: int = 100, max_age_s: float = 30 * 86400) -> dict[str, dict | None]:
    """Polymer entities per entry: UniProt IDs, organisms, polymer type, Pfam/membrane annotations.

    Cached as one JSON file (large; not committed). Unknown entries map to None.
    """
    path = CACHE_DIR / "entry_entities.json"
    cache: dict[str, dict | None] = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    todo = sorted({p.upper() for p in pdb_ids} - cache.keys())
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        q = "query($ids: [String!]!) { entries(entry_ids: $ids) { %s } }" % _ENTITY_FIELDS
        data = _post_json(GRAPHQL, {"query": q, "variables": {"ids": chunk}})
        found = {e["rcsb_id"]: e for e in (data.get("data") or {}).get("entries") or [] if e}
        for pid in chunk:
            cache[pid] = found.get(pid)
        if (i // batch) % 10 == 0 or i + batch >= len(todo):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    return {p: cache.get(p.upper()) for p in pdb_ids}


def entry_keywords(pdb_ids: list[str], batch: int = 200, max_age_s: float = 30 * 86400) -> dict[str, list[str] | None]:
    """UniProt keyword IDs (e.g. KW-0812 Transmembrane) over all protein entities of each entry."""
    path = CACHE_DIR / "entry_keywords.json"
    cache: dict[str, list[str] | None] = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    todo = sorted({p.upper() for p in pdb_ids} - cache.keys())
    q = ("query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id polymer_entities { "
         "uniprots { rcsb_uniprot_keyword { id } } } } }")
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        data = _post_json(GRAPHQL, {"query": q, "variables": {"ids": chunk}})
        found = {}
        for e in (data.get("data") or {}).get("entries") or []:
            if not e:
                continue
            kws = {k["id"] for pe in e.get("polymer_entities") or [] for u in pe.get("uniprots") or []
                   for k in u.get("rcsb_uniprot_keyword") or []}
            found[e["rcsb_id"]] = sorted(kws)
        for pid in chunk:
            cache[pid] = found.get(pid)
        if (i // batch) % 10 == 0 or i + batch >= len(todo):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    return {p: cache.get(p.upper()) for p in pdb_ids}


def entry_extras(pdb_ids: list[str], batch: int = 150, max_age_s: float = 30 * 86400) -> dict[str, dict | None]:
    """Per entry: UniProt keywords per accession, and when its citation was added.

    {"kw": {uniprot: [KW-...]}, "released": date, "db_ref_updates": [dates of
    'Database references' revisions after release]} -- the PDB adds or updates
    the primary citation in such revisions, so a citation that arrived later
    was not known in the structure's release week.
    """
    path = CACHE_DIR / "entry_extras.json"
    cache: dict[str, dict | None] = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    todo = sorted({p.upper() for p in pdb_ids} - cache.keys())
    q = ("query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id "
         "rcsb_accession_info { initial_release_date } "
         "pdbx_audit_revision_history { ordinal revision_date } "
         "pdbx_audit_revision_group { revision_ordinal group } "
         "polymer_entities { uniprots { rcsb_id rcsb_uniprot_keyword { id } } } } }")
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        data = _post_json(GRAPHQL, {"query": q, "variables": {"ids": chunk}})
        found = {}
        for e in (data.get("data") or {}).get("entries") or []:
            if not e:
                continue
            rel = ((e.get("rcsb_accession_info") or {}).get("initial_release_date") or "")[:10]
            dates = {h["ordinal"]: (h.get("revision_date") or "")[:10] for h in e.get("pdbx_audit_revision_history") or []}
            refs = sorted({dates.get(g["revision_ordinal"], "") for g in e.get("pdbx_audit_revision_group") or []
                           if g.get("group") == "Database references"} - {""})
            kw = {}
            for pe in e.get("polymer_entities") or []:
                for u in pe.get("uniprots") or []:
                    kw[u["rcsb_id"]] = sorted({k["id"] for k in u.get("rcsb_uniprot_keyword") or []})
            found[e["rcsb_id"]] = {"released": rel, "db_ref_updates": [d for d in refs if d > rel], "kw": kw}
        for pid in chunk:
            cache[pid] = found.get(pid)
        if (i // batch) % 10 == 0 or i + batch >= len(todo):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    return {p: cache.get(p.upper()) for p in pdb_ids}


def entry_families(pdb_ids: list[str], batch: int = 150, max_age_s: float = 30 * 86400) -> dict[str, dict | None]:
    """Depositing authors and UniProt's InterPro families per accession.

    Both exist on release day (unlike the PDB's own Pfam/InterPro labels, which
    arrive later), so new and old structures can be compared the same way.
    {"authors": [...], "ipr": {uniprot: [IPR...]}}
    """
    path = CACHE_DIR / "entry_families.json"
    cache: dict[str, dict | None] = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    todo = sorted({p.upper() for p in pdb_ids} - cache.keys())
    q = ("query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id audit_author { name } "
         "polymer_entities { uniprots { rcsb_id rcsb_uniprot_annotation { type annotation_id } } } } }")
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        data = _post_json(GRAPHQL, {"query": q, "variables": {"ids": chunk}})
        found = {}
        for e in (data.get("data") or {}).get("entries") or []:
            if not e:
                continue
            ipr = {}
            for pe in e.get("polymer_entities") or []:
                for u in pe.get("uniprots") or []:
                    ipr[u["rcsb_id"]] = sorted({a["annotation_id"] for a in u.get("rcsb_uniprot_annotation") or []
                                                if a.get("type") == "InterPro"})
            found[e["rcsb_id"]] = {"authors": [a["name"] for a in e.get("audit_author") or [] if a.get("name")], "ipr": ipr}
        for pid in chunk:
            cache[pid] = found.get(pid)
        if (i // batch) % 10 == 0 or i + batch >= len(todo):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    return {p: cache.get(p.upper()) for p in pdb_ids}


def entry_details(pdb_ids: list[str], batch: int = 100, max_age_s: float = 7 * 86400) -> dict[str, dict | None]:
    """What an mpstruc-style listing needs: entity names, organisms with lineage,
    expression hosts, bound ligands and the full primary citation."""
    path = CACHE_DIR / "entry_details.json"
    cache: dict[str, dict | None] = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    todo = sorted({p.upper() for p in pdb_ids} - cache.keys())
    q = ("query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id struct { title } exptl { method } "
         "rcsb_entry_info { resolution_combined } "
         "rcsb_primary_citation { rcsb_authors title journal_abbrev journal_volume journal_issue page_first page_last "
         "year pdbx_database_id_DOI pdbx_database_id_PubMed } "
         "polymer_entities { rcsb_polymer_entity { pdbx_description } "
         "entity_poly { rcsb_entity_polymer_type } "
         "rcsb_polymer_entity_container_identifiers { uniprot_ids } "
         "rcsb_entity_source_organism { ncbi_scientific_name taxonomy_lineage { name } } "
         "rcsb_entity_host_organism { ncbi_scientific_name } } "
         "nonpolymer_entities { nonpolymer_comp { chem_comp { id name } } } } }")
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        data = _post_json(GRAPHQL, {"query": q, "variables": {"ids": chunk}})
        found = {e["rcsb_id"]: e for e in (data.get("data") or {}).get("entries") or [] if e}
        for pid in chunk:
            cache[pid] = found.get(pid)
        if (i // batch) % 10 == 0 or i + batch >= len(todo):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    return {p: cache.get(p.upper()) for p in pdb_ids}
