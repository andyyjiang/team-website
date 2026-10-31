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
