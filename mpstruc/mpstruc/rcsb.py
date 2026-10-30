"""Minimal cached RCSB PDB client (holdings only, for now)."""
from __future__ import annotations

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
        except urllib.error.URLError:
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
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
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
