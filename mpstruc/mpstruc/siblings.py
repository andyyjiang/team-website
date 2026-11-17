"""Same-paper structures released before the weekly window.

The weekly run only looks at PDB releases of the last 26 weeks. A paper's
structures are often released months apart (9WFO, GGCX with proP1, came out
on 15 July 2026; the six other GGCX complexes of the same paper on 1 and 22
April), and mpstruc lists a paper's structures together: one entry with the
others as related lines. So for every paper of a candidate the run asks the
PDB for *all* entries that cite it, and takes in those that are not in the
window, not in mpstruc and not decided by a curator: the "siblings". They are
scored and placed with the same code as every other candidate, so the study
groups together under one lead (the lowest PDB ID, see automation.batch_leads),
and each carries

    "sibling": {"of": the paper's DOI (or "pmid:<n>"), "with": [the paper's
                structures inside the window], "reason": a sentence}

In official mpstruc, 50 of the 2,337 entries that list same-paper related
structures span releases more than 26 weeks apart, so curators do group them.

The functions that decide (paper_key, papers, select, keep, note) are pure;
find_entries talks to the PDB search API (batched: one request per BATCH
papers, cached in cache/paper_entries.json for MAX_AGE_S) and tolerates
failures: a paper that could not be looked up simply has no siblings this run.
"""
from __future__ import annotations

import datetime as dt
import json
import time
import urllib.request
from pathlib import Path
from typing import Callable

from . import awaiting, rcsb

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
CACHE_NAME = "paper_entries.json"
MAX_AGE_S = 7 * 86400  # the weekly run asks again; a second run in the same week does not
BATCH = 50  # papers per search request
PAUSE_S = 0.25  # between requests
DOI_ATTR = "rcsb_primary_citation.pdbx_database_id_DOI"
# The search API does not index rcsb_primary_citation.pdbx_database_id_PubMed; this one it does.
PMID_ATTR = "rcsb_pubmed_container_identifiers.pubmed_id"
STAT = "same-paper structures from before the window"


def paper_key(doi: str | None, pmid=None) -> str | None:
    """'doi:10.1038/…' (lower case), or 'pmid:123' when the PDB gives no DOI, or None."""
    doi = str(doi or "").strip().lower()
    if doi:
        return f"doi:{doi}"
    pm = str(pmid or "").strip()
    return f"pmid:{pm}" if pm.isdigit() and int(pm) > 0 else None


def citation_key(cit: dict | None) -> str | None:
    cit = cit or {}
    return paper_key(cit.get("pdbx_database_id_DOI"), cit.get("pdbx_database_id_PubMed"))


def papers(cands: list[dict]) -> dict[str, list[str]]:
    """{paper: [its structures among the candidates]} for candidates whose paper is out.

    cands: [{"pdb", "doi", "pmid"[, "awaiting", "sibling"]}]. Structures waiting for
    their paper (mpstruc does not list them yet) and siblings found by an earlier
    call are left out: only a structure of the window brings its paper in."""
    out: dict[str, list[str]] = {}
    for c in cands:
        if c.get("awaiting") or c.get("sibling"):
            continue
        key = paper_key(c.get("doi"), c.get("pmid"))
        if key:
            out.setdefault(key, []).append(c["pdb"].upper())
    return {k: sorted(set(v)) for k, v in sorted(out.items())}


def search_query(keys: list[str]) -> dict:
    """One search request for every entry whose primary citation is one of these papers."""
    dois = [k[4:] for k in keys if k.startswith("doi:")]
    pmids = [int(k[5:]) for k in keys if k.startswith("pmid:")]
    nodes = []
    if dois:
        nodes.append({"type": "terminal", "service": "text",
                      "parameters": {"attribute": DOI_ATTR, "operator": "in", "value": dois}})
    if pmids:
        nodes.append({"type": "terminal", "service": "text",
                      "parameters": {"attribute": PMID_ATTR, "operator": "in", "value": pmids}})
    query = nodes[0] if len(nodes) == 1 else {"type": "group", "logical_operator": "or", "nodes": nodes}
    return {"query": query, "return_type": "entry", "request_options": {"return_all_hits": True}}


def _search(query: dict) -> list[str]:
    req = urllib.request.Request(SEARCH, data=json.dumps(query).encode(),
                                 headers={"Content-Type": "application/json", "User-Agent": "mpstruc-curation-tool"})
    with urllib.request.urlopen(req, timeout=60) as r:
        body = r.read()
    return sorted(x["identifier"] for x in json.loads(body).get("result_set", [])) if body else []  # 204: no hits


def find_entries(keys: list[str], search: Callable[[dict], list[str]] | None = None,
                 citations: Callable[[list[str]], dict[str, dict | None]] | None = None,
                 cache_path: Path | None = None, now: float | None = None, pause: float = PAUSE_S,
                 max_age_s: float = MAX_AGE_S) -> tuple[dict[str, list[str]], dict[str, dict | None]]:
    """Every PDB entry citing each paper: ({paper: [pdb, ...]}, {pdb: its citation}).

    Papers looked up less than max_age_s ago come from the cache. The others are
    asked BATCH at a time; the entries found are then matched to their paper by
    their own citation (one batched GraphQL query, awaiting.fresh_citations). A
    batch that fails is left out of the result (and of the cache), so a later run
    asks again; the citations returned are those fetched in this call only."""
    search = search or _search
    citations = citations or awaiting.fresh_citations
    cache_path = cache_path or rcsb.CACHE_DIR / CACHE_NAME
    now = time.time() if now is None else now
    try:
        cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    except ValueError:
        cache = {}
    keys = sorted(set(keys))
    out = {k: list(cache[k]["ids"]) for k in keys if k in cache and now - cache[k].get("at", 0) < max_age_s}
    todo = [k for k in keys if k not in out]
    cites: dict[str, dict | None] = {}
    failed = 0
    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        if i and pause:
            time.sleep(pause)
        try:
            ids = search(search_query(chunk))
            got = citations(ids) if ids else {}
        except Exception as err:  # the PDB is unreachable or refuses: these papers get no siblings this run
            failed += len(chunk)
            print(f"siblings: could not look up {len(chunk)} papers ({err})")
            continue
        found: dict[str, list[str]] = {k: [] for k in chunk}
        for pdb, cit in got.items():
            cit = cit or {}
            # A paper is found by its DOI, or by its PubMed id when a candidate has no DOI.
            for key in {paper_key(cit.get("pdbx_database_id_DOI")), paper_key(None, cit.get("pdbx_database_id_PubMed"))}:
                if key in found:
                    found[key].append(pdb.upper())
        cites.update({p.upper(): c for p, c in got.items()})
        for k, v in found.items():
            out[k] = sorted(set(v))
            cache[k] = {"at": now, "ids": out[k]}
    if len(todo) > failed:
        # Papers nobody asked about for a long time are dropped, so the file stays small.
        cache = {k: v for k, v in cache.items() if now - v.get("at", 0) < 8 * max_age_s}
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache))
    return out, cites


def select(by_paper: dict[str, list[str]], found: dict[str, list[str]], released: dict[str, str], since: str,
           have: set[str], listed: set[str], decided: set[str],
           unpublished: Callable[[str], bool] = lambda pdb: False) -> dict[str, dict]:
    """The siblings: {pdb: {"paper": key, "with": [the paper's candidates in the window]}}.

    by_paper: papers(); found: find_entries(); released: {pdb: release date}; since:
    first day of the window. Left out: structures the run already has (have: the
    window's releases and those tracked while they wait for a paper), those in
    official mpstruc (listed: an entry, a related code or named in a description),
    those a curator decided (decided: autoreject.decided_ids, which leaves out 'Not
    published yet' rejections, so those follow the usual rules), those not released
    before the window (or with no release date), and those whose own PDB record does
    not show the paper as out (unpublished)."""
    out: dict[str, dict] = {}
    for key, mates in by_paper.items():
        for pdb in found.get(key) or []:
            pdb = pdb.upper()
            rel = released.get(pdb) or ""
            if pdb in have or pdb in listed or pdb in decided or not rel or rel >= since or unpublished(pdb):
                continue
            if pdb not in out:  # (a paper found by DOI and by PubMed id is the same paper)
                out[pdb] = {"paper": key, "with": list(mates)}
    return dict(sorted(out.items()))


def keep(sibs: dict[str, dict], scores: dict[str, float], floor: float, cores: dict[str, set],
         kept: set[str]) -> dict[str, str]:
    """Which siblings become candidates: {pdb: why}.

    The same floor as every other structure (scores: membrane-protein probability),
    except that a sibling with exactly the proteins of a candidate of the same paper
    that was kept is kept too: the same protein in the same paper belongs to the
    same mpstruc entry whatever its own score says. cores: {pdb: membrane proteins}
    of the siblings and of those candidates (kept)."""
    out = {}
    for pdb, s in sibs.items():
        if pdb not in scores:
            continue  # the PDB gave nothing to score
        if scores[pdb] >= floor:
            out[pdb] = "score"
            continue
        core = cores.get(pdb) or set()
        same = [m for m in s["with"] if m in kept and core and cores.get(m) == core]
        if same:
            out[pdb] = f"same protein as {same[0]}"
    return out


def weeks_text(since: str, until: str) -> str:
    """'26-week' for the usual window."""
    try:
        days = (dt.date.fromisoformat(until) - dt.date.fromisoformat(since)).days
    except ValueError:
        return "release"
    return f"{round(days / 7)}-week"


def note(sib: dict, window: list[str]) -> dict:
    """The candidate's "sibling" field."""
    mates = sib["with"]
    names = mates[0] if len(mates) == 1 else f"{mates[0]} and {len(mates) - 1} more"
    key = sib["paper"]
    return {"of": key[4:] if key.startswith("doi:") else key, "with": list(mates),
            "reason": f"released before the {weeks_text(*window)} window; same paper as {names}"}


def stale(cached: dict[str, dict | None], cites: dict[str, dict | None]) -> list[str]:
    """Entries whose cached PDB record names another paper (or none) than the citation just fetched."""
    return sorted(p for p, cit in cites.items()
                  if p in cached and citation_key((cached[p] or {}).get("rcsb_primary_citation")) != citation_key(cit))


def discover(cands: list[dict], window: list[str], have: set[str], listed: set[str], decided: set[str],
             find=find_entries) -> tuple[dict[str, dict], dict]:
    """select() for these candidates, with the PDB lookups: (siblings, counts).

    Each sibling also gets "released". A sibling whose own PDB record still says
    'To be published' or names only a preprint is left for the week that changes.
    Entries whose citation changed since they were cached are dropped from the
    local PDB caches (awaiting.forget), so they are placed by today's paper."""
    by_paper = papers(cands)
    t0 = time.time()
    found, cites = find(sorted(by_paper))
    others = sorted({p for k in by_paper for p in found.get(k) or []} - have - listed - decided)
    released = {}
    if others:
        missing = [p for p in others if p not in cites]  # (found in the cache: their citation was not fetched now)
        if missing:
            try:
                cites.update(awaiting.fresh_citations(missing))
            except Exception as err:
                print(f"siblings: could not check {len(missing)} citations ({err})")
        try:
            # An old structure may have been cached before the PDB linked its paper: it would
            # then not be placed with the paper's other structures.
            awaiting.forget(stale(rcsb.entry_entities(others), cites))
        except Exception as err:
            print(f"siblings: could not refresh the PDB caches ({err})")
        ents = rcsb.entry_entities(others)
        released = {p: (((ents.get(p) or {}).get("rcsb_accession_info") or {}).get("initial_release_date") or "")[:10]
                    for p in others}
    sibs = select(by_paper, found, released, window[0], have, listed, decided,
                  unpublished=lambda p: p not in cites or bool(awaiting.citation_status(cites[p])))
    for p, s in sibs.items():
        s["released"] = released[p]
    counts = {"papers": len(by_paper), "papersLookedUp": sum(k in found for k in by_paper),
              "papersWithSiblings": len({s["paper"] for s in sibs.values()}), "found": len(sibs),
              "seconds": round(time.time() - t0, 1)}
    return sibs, counts
