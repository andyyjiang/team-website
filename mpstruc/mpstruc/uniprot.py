"""UniProt membrane segments (transmembrane and intramembrane regions)."""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request

from .rcsb import CACHE_DIR

SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def tm_segments(accessions: list[str], batch: int = 100, max_age_s: float = 30 * 86400) -> dict[str, list[list[int]]]:
    """{accession: [[first, last], ...]} membrane-crossing segments in UniProt numbering ([] if none)."""
    path = CACHE_DIR / "uniprot_tm.json"
    cache: dict[str, list] = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    todo = sorted({a for a in accessions if a and re.fullmatch(r"[A-Z0-9]{6,10}", a)} - cache.keys())
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        query = " OR ".join(f"accession:{a}" for a in chunk)
        url = SEARCH + "?" + urllib.parse.urlencode({"query": query, "fields": "accession,ft_transmem,ft_intramem",
                                                     "format": "tsv", "size": 500})
        for attempt in range(5):
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "mpstruc-curation-tool"}),
                                            timeout=120) as r:
                    text = r.read().decode()
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(2 ** (attempt + 1))
        for a in chunk:
            cache[a] = []
        for line in text.splitlines()[1:]:
            acc, *cols = line.split("\t")
            segs = [[int(b), int(e)] for col in cols for b, e in re.findall(r"(?:TRANSMEM|INTRAMEM) (\d+)\.\.(\d+)", col)]
            cache[acc] = segs
        if (i // batch) % 20 == 0 or i + batch >= len(todo):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    return {a: cache.get(a, []) for a in accessions}
