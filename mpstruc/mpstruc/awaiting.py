"""Structures waiting for their paper.

mpstruc lists a structure only once it is released *and* its peer-reviewed
paper is out (with page numbers, in PubMed). Many new PDB entries still cite
"To be published", or only a preprint. The weekly run marks those
suggestions with

    "awaiting": {"since": first run date that saw it waiting, "reason": "To be published" | "preprint only"}

so the site keeps them out of the default "To decide" list, and remembers
them in data/awaiting.json:

    {"updated": run date,
     "awaiting":  {pdb: {"since": date, "released": date, "reason": ...}},
     "published": {pdb: {"on": run date that first saw the paper, "released": date}}}

Every ID in that file is re-checked by later runs even after it falls out of
the 26-week release window, so a structure comes back as a normal suggestion
(with "publishedOn") whenever its paper appears, however late. The
'published' map keeps that date for about 26 weeks.

The functions that decide (citation_status, update) are pure; the network
and cache helpers below them are used only by automation.main().
"""
from __future__ import annotations

import datetime as dt
import json
import re
from pathlib import Path

KEEP_PUBLISHED_DAYS = 182  # ~26 weeks
# A structure still without a paper this long after its release (or, if the
# release date is unknown, after we first saw it waiting) is given up on:
# many PDB entries never get a paper, and re-checking them forever would grow
# the weekly run and the curators' "Awaiting publication" list without end.
GIVE_UP_DAYS = 730  # ~2 years

TO_BE_PUBLISHED = "To be published"
PREPRINT_ONLY = "preprint only"

# Journal names the PDB uses for "no paper yet" (seen in any capitalisation).
_UNPUBLISHED = {"to be published", "not published", "unpublished", "in preparation", "submitted",
                "to be submitted", "manuscript in preparation"}
# Preprint servers, matched against the journal name with spaces and dots removed.
_PREPRINT = ("biorxiv", "medrxiv", "chemrxiv", "arxiv", "researchsquare", "ressq", "ssrn",
             "preprints", "preprint")
# Preprint DOIs: bioRxiv/medRxiv (10.1101/2024.01.02.123456, not Cold Spring
# Harbor journals, which share the prefix; also bioRxiv's newer 10.64898/),
# Research Square, arXiv, SSRN, ChemRxiv, Preprints.org.
_PREPRINT_DOI = re.compile(r"^10\.(1101/\d{4}\.\d{2}\.\d{2}\.|64898/|21203/|48550/|2139/|26434/|20944/)", re.I)


def citation_status(cit: dict | None) -> str | None:
    """Why a structure is not ready for mpstruc yet, or None if its paper is out.

    cit is the PDB's rcsb_primary_citation (None when the entry cites nothing).
    """
    if not cit:
        return TO_BE_PUBLISHED
    journal = " ".join(str(cit.get("journal_abbrev") or "").split()).strip(" .").lower()
    if journal in _UNPUBLISHED:
        return TO_BE_PUBLISHED
    squashed = re.sub(r"[\s.]", "", journal)
    if any(squashed.startswith(p) for p in _PREPRINT):
        return PREPRINT_ONLY
    doi = str(cit.get("pdbx_database_id_DOI") or "").strip()
    if not journal:
        if doi and _PREPRINT_DOI.match(doi):
            return PREPRINT_ONLY
        if not doi and not cit.get("pdbx_database_id_PubMed"):
            return TO_BE_PUBLISHED
    return None


def empty() -> dict:
    return {"updated": None, "awaiting": {}, "published": {}}


def load(path: Path) -> dict:
    """The saved record, or an empty one. Accepts a plain {pdb: {...}} map too."""
    try:
        raw = json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return empty()
    if not isinstance(raw, dict):
        return empty()
    if "awaiting" not in raw and "published" not in raw:
        raw = {"awaiting": {k: v for k, v in raw.items() if isinstance(v, dict)}}
    rec = empty()
    rec["updated"] = raw.get("updated")
    rec["awaiting"] = {k.upper(): dict(v) for k, v in (raw.get("awaiting") or {}).items() if isinstance(v, dict)}
    rec["published"] = {k.upper(): dict(v) for k, v in (raw.get("published") or {}).items() if isinstance(v, dict)}
    return rec


def save(path: Path, rec: dict) -> None:
    out = {"updated": rec.get("updated"),
           "awaiting": dict(sorted(rec["awaiting"].items())),
           "published": dict(sorted(rec["published"].items()))}
    Path(path).write_text(json.dumps(out, indent=1) + "\n")


def tracked(rec: dict) -> set[str]:
    """IDs every run must re-check, whatever their release date."""
    return set(rec["awaiting"]) | set(rec["published"])


def drop_stale(rec: dict, today: str) -> tuple[dict, list[str]]:
    """Stop tracking structures that have waited longer than GIVE_UP_DAYS.

    Returns the record without them and their IDs. They are no longer
    re-checked; if one is still inside the release window it simply stays
    marked as waiting for this run and is forgotten after that.
    """
    cutoff = (dt.date.fromisoformat(today) - dt.timedelta(days=GIVE_UP_DAYS)).isoformat()
    keep, gone = {}, []
    for pdb, v in (rec.get("awaiting") or {}).items():
        start = str(v.get("released") or v.get("since") or "")[:10]
        if start and start < cutoff:
            gone.append(pdb)
        else:
            keep[pdb] = v
    return {**rec, "awaiting": keep}, sorted(gone)


def with_tracked(ids: list[str], rec: dict) -> list[str]:
    """The window's releases plus tracked structures released before it."""
    have = {p.upper() for p in ids}
    return list(ids) + sorted(p for p in tracked(rec) if p not in have)


def update(cands: list[dict], citations: dict[str, dict | None | bool], rec: dict, today: str) -> tuple[dict, dict]:
    """Mark candidates waiting for their paper and return the new record and counts.

    citations: pdb -> rcsb_primary_citation (dict or None for "cites nothing");
    a pdb missing from it (PDB data unavailable) keeps its previous state.
    Sets c["awaiting"] or c["publishedOn"] on the candidates in place.
    Tracked structures that are no longer candidates (approved, rejected,
    added to mpstruc, or obsoleted) are dropped from the record.
    """
    old_wait, old_pub = rec.get("awaiting") or {}, rec.get("published") or {}
    cutoff = (dt.date.fromisoformat(today) - dt.timedelta(days=KEEP_PUBLISHED_DAYS)).isoformat()
    wait, pub = {}, {}
    newly = 0
    for c in cands:
        pdb = c["pdb"].upper()
        c.pop("awaiting", None)
        c.pop("publishedOn", None)
        if pdb in citations:
            reason = citation_status(citations[pdb])
        elif pdb in old_wait:
            reason = old_wait[pdb].get("reason") or TO_BE_PUBLISHED
        else:
            reason = None
        if reason:
            since = (old_wait.get(pdb) or {}).get("since") or today
            wait[pdb] = {"since": since, "released": c.get("released"), "reason": reason}
            c["awaiting"] = {"since": since, "reason": reason}
            continue
        on = (old_pub.get(pdb) or {}).get("on")
        if not on and pdb in old_wait:
            on = today
            newly += 1
        if on and on >= cutoff:
            pub[pdb] = {"on": on, "released": c.get("released")}
            c["publishedOn"] = on
    new = {"updated": today, "awaiting": wait, "published": pub}
    stats = {"awaiting": len(wait), "newly_published": newly}
    return new, stats


# ------------------------------------------------------------- PDB helpers

# Cache files that hold per-entry PDB data, with the age (seconds) after which
# rcsb.py stops trusting them (it goes by file modification time).
_CACHES = {"entry_details.json": 7 * 86400, "entry_entities.json": 30 * 86400,
           "entry_summaries.json": 7 * 86400, "entry_extras.json": 30 * 86400}


def fresh_citations(ids: list[str], batch: int = 200) -> dict[str, dict | None]:
    """Primary citation per entry straight from the PDB (no cache), batched.

    Entries the PDB does not return are left out, so their state is kept.
    """
    from . import rcsb
    out: dict[str, dict | None] = {}
    ids = sorted({p.upper() for p in ids})
    q = ("query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id "
         "rcsb_primary_citation { title journal_abbrev journal_volume page_first year "
         "pdbx_database_id_DOI pdbx_database_id_PubMed } } }")
    for i in range(0, len(ids), batch):
        data = rcsb._post_json(rcsb.GRAPHQL, {"query": q, "variables": {"ids": ids[i:i + batch]}})
        for e in (data.get("data") or {}).get("entries") or []:
            if e:
                out[e["rcsb_id"]] = e.get("rcsb_primary_citation")
    return out


def _cite_key(cit: dict | None) -> tuple:
    c = cit or {}
    return (c.get("journal_abbrev"), c.get("pdbx_database_id_DOI"), c.get("pdbx_database_id_PubMed"), c.get("page_first"))


def forget_changed(fresh: dict[str, dict | None]) -> list[str]:
    """Drop entries whose citation changed from the local PDB caches, so the
    run fetches their new paper details (and the draft listing cites it)."""
    import os
    import time
    from . import rcsb

    def expired(f: Path) -> bool:
        return time.time() - f.stat().st_mtime >= _CACHES[f.name]

    path = rcsb.CACHE_DIR / "entry_details.json"
    try:
        # An expired cache is refetched in full anyway, so nothing to forget.
        cached = json.loads(path.read_text()) if path.exists() and not expired(path) else {}
    except ValueError:
        cached = {}
    changed = [p for p, cit in fresh.items()
               if p in cached and _cite_key((cached.get(p) or {}).get("rcsb_primary_citation")) != _cite_key(cit)]
    if not changed:
        return []
    for name in _CACHES:
        f = rcsb.CACHE_DIR / name
        if not f.exists() or expired(f):
            continue  # rewriting an expired file would make it look new again
        try:
            cache = json.loads(f.read_text())
        except ValueError:
            continue
        if any(p in cache for p in changed):
            for p in changed:
                cache.pop(p, None)
            st = f.stat()
            f.write_text(json.dumps(cache))
            # Keep the original age, so pruning never extends the cache's life.
            os.utime(f, ns=(st.st_atime_ns, st.st_mtime_ns))
    return changed


def citations_for(cands: list[dict], fresh: dict[str, dict | None]) -> dict[str, dict | None]:
    """Citation per candidate: fresh where already fetched, else the cached
    entry details; anything that looks unpublished in the cache is checked
    fresh (one batched query) so an old cache never hides a new paper.

    Entries whose re-checked citation differs from the cache are dropped from
    the local caches (forget_changed), so call this before drafting listings:
    the drafts then cite the new paper rather than "To be published"."""
    from . import rcsb
    ids = [c["pdb"] for c in cands]
    det = rcsb.entry_details(ids)
    out = {p: fresh[p] for p in ids if p in fresh}
    recheck = []
    for p in ids:
        if p in out or det.get(p) is None:
            continue
        cit = det[p].get("rcsb_primary_citation")
        out[p] = cit
        if citation_status(cit):
            recheck.append(p)
    if recheck:
        try:
            got = fresh_citations(recheck)
        except Exception as err:  # network trouble: the cached answer stands
            print(f"awaiting: could not re-check {len(recheck)} citations ({err}); using cached ones")
            got = {}
        out.update(got)
        forget_changed(got)
    return out
