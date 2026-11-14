"""Minimal cached PubMed (NCBI E-utilities esummary) client."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request

from .rcsb import CACHE_DIR

ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def summaries(pmids: list[str], batch: int = 200, max_age_s: float = 30 * 86400) -> dict[str, dict | None]:
    """{pmid: {"title", "source", "year", "doi"}} for valid numeric PMIDs; None if PubMed has no record."""
    path = CACHE_DIR / "pubmed_summaries.json"
    cache: dict[str, dict | None] = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    want = sorted({p for p in pmids if p and p.isdigit()})
    todo = [p for p in want if p not in cache]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        body = urllib.parse.urlencode({"db": "pubmed", "retmode": "json", "id": ",".join(chunk)}).encode()
        req = urllib.request.Request(ESUMMARY, data=body, headers={"User-Agent": "mpstruc-curation-tool"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    res = json.load(r).get("result", {})
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2 ** (attempt + 1))
        for p in chunk:
            x = res.get(p)
            if not x or x.get("error"):
                cache[p] = None
                continue
            doi = next((a["value"] for a in x.get("articleids", []) if a.get("idtype") == "doi"), None)
            cache[p] = {"title": x.get("title"), "source": x.get("source"),
                        "year": (x.get("pubdate") or "")[:4], "doi": doi}
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cache))
        time.sleep(0.4)  # NCBI allows 3 requests/s without an API key
    return {p: cache.get(p) for p in pmids}


ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CROSSREF = "https://api.crossref.org/works"
MISS_RETRY_S = 6 * 86400  # a paper not found yet is looked for again by next week's run


def _norm_title(t: str | None) -> str:
    import html
    import re
    t = html.unescape(re.sub(r"<[^>]+>", "", t or ""))
    return re.sub(r"[^a-z0-9]+", " ", t.lower()).strip()


def _fetch(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "mpstruc-curation-tool"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 ** (attempt + 1))
    return {}


def find_paper(title: str | None, year: str | int | None) -> dict | None:
    """A published paper with exactly this title, from PubMed or else Crossref:
    {"pmid": str | None, "doi": str | None, "source": "pubmed" | "crossref"}, or None.

    For PDB citations that name a journal but carry no DOI or PubMed ID: the
    PDB fills those in late, so a journal name alone does not say whether the
    paper is out. The title must match exactly (ignoring case and punctuation)
    and the year must be within one of the PDB's, so a placeholder title
    ('Membrane Protein Transporter') never matches an older paper. Results are
    cached; a miss is looked for again after MISS_RETRY_S."""
    norm = _norm_title(title)
    if len(norm.split()) < 4:
        return None  # too short to identify a paper
    try:
        y = int(str(year)[:4])
    except (TypeError, ValueError):
        y = None
    key = f"{norm}|{y}"
    path = CACHE_DIR / "paper_lookup.json"
    cache = json.loads(path.read_text()) if path.exists() else {}
    hit = cache.get(key)
    if hit is not None and (hit.get("found") or time.time() - hit.get("at", 0) < MISS_RETRY_S):
        return hit.get("found")

    def year_ok(got) -> bool:
        try:
            return y is None or abs(int(str(got)[:4]) - y) <= 1
        except (TypeError, ValueError):
            return False

    found = None
    q = urllib.parse.urlencode({"db": "pubmed", "retmode": "json", "retmax": 5, "term": f"{title}[ti]"})
    ids = (_fetch(f"{ESEARCH}?{q}").get("esearchresult") or {}).get("idlist") or []
    time.sleep(0.4)
    if ids:
        for pmid, s in summaries(ids, max_age_s=0).items():
            if s and _norm_title(s.get("title")) == norm and year_ok(s.get("year")):
                found = {"pmid": pmid, "doi": s.get("doi"), "source": "pubmed"}
                break
    if not found:
        q = urllib.parse.urlencode({"query.bibliographic": title, "rows": 5, "filter": "type:journal-article",
                                    "select": "DOI,title,published"})
        for it in (_fetch(f"{CROSSREF}?{q}").get("message") or {}).get("items") or []:
            got = ((it.get("published") or {}).get("date-parts") or [[None]])[0][0]
            if _norm_title((it.get("title") or [""])[0]) == norm and year_ok(got):
                found = {"pmid": None, "doi": it.get("DOI"), "source": "crossref"}
                break
    cache[key] = {"found": found, "at": time.time()}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache))
    return found
