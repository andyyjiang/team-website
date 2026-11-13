"""UniProt membrane segments and named domains, batched from rest.uniprot.org.

- tm_segments: transmembrane and intramembrane regions (ft_transmem, ft_intramem).
- membrane_features: the same segments plus UniProt's Domain and Region
  features with their names (ft_domain, ft_region). These are the bars in the
  feature viewer at the bottom of an RCSB PDB page. Outer-membrane beta-barrels
  often have no transmembrane features but do have a domain such as
  'TBDR beta-barrel' or 'Autotransporter' (see membrane.barrel_domains).
"""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request

from .rcsb import CACHE_DIR

SEARCH = "https://rest.uniprot.org/uniprotkb/search"
_SEG = re.compile(r"(?:TRANSMEM|INTRAMEM) (\d+)\.\.(\d+)")
_NAMED = re.compile(r'(?:DOMAIN|REGION) (\d+)\.\.(\d+); /note="([^"]*)"')


def _query(fields: str, chunk: list[str]) -> str:
    query = " OR ".join(f"accession:{a}" for a in chunk)
    url = SEARCH + "?" + urllib.parse.urlencode({"query": query, "fields": fields, "format": "tsv", "size": 500})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "mpstruc-curation-tool"}),
                                        timeout=120) as r:
                return r.read().decode()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** (attempt + 1))
    raise RuntimeError("unreachable")


def _batched(name: str, fields: str, parse, accessions: list[str], batch: int, max_age_s: float, empty) -> dict:
    path = CACHE_DIR / name
    cache: dict = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    todo = sorted({a for a in accessions if a and re.fullmatch(r"[A-Z0-9]{6,10}", a)} - cache.keys())
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        text = _query(fields, chunk)
        for a in chunk:
            cache[a] = empty()
        for line in text.splitlines()[1:]:
            acc, *cols = line.split("\t")
            cache[acc] = parse(cols)
        if (i // batch) % 20 == 0 or i + batch >= len(todo):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    return {a: cache.get(a, empty()) for a in accessions}


def tm_segments(accessions: list[str], batch: int = 100, max_age_s: float = 30 * 86400) -> dict[str, list[list[int]]]:
    """{accession: [[first, last], ...]} membrane-crossing segments in UniProt numbering ([] if none)."""
    return _batched("uniprot_tm.json", "accession,ft_transmem,ft_intramem",
                    lambda cols: [[int(b), int(e)] for col in cols for b, e in _SEG.findall(col)],
                    accessions, batch, max_age_s, list)


def membrane_features(accessions: list[str], batch: int = 100, max_age_s: float = 30 * 86400) -> dict[str, dict]:
    """{accession: {"tm": [[first, last], ...], "named": [[first, last, name], ...]}}.

    "tm": transmembrane/intramembrane segments; "named": Domain and Region features."""
    def parse(cols):
        tm, dom_reg = cols[:2], cols[2:]
        return {"tm": [[int(b), int(e)] for col in tm for b, e in _SEG.findall(col)],
                "named": [[int(b), int(e), n] for col in dom_reg for b, e, n in _NAMED.findall(col)]}
    return _batched("uniprot_membrane.json", "accession,ft_transmem,ft_intramem,ft_domain,ft_region",
                    parse, accessions, batch, max_age_s, lambda: {"tm": [], "named": []})


ACCESSIONS = "https://rest.uniprot.org/uniprotkb/accessions"


def _get(url: str) -> bytes:
    for attempt in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "mpstruc-curation-tool"}),
                                        timeout=120) as r:
                return r.read()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** (attempt + 1))
    raise RuntimeError("unreachable")


def _protein_name(desc: dict | None) -> str:
    desc = desc or {}
    for names in ([desc.get("recommendedName")], desc.get("submissionNames") or [], desc.get("alternativeNames") or []):
        for n in names:
            value = ((n or {}).get("fullName") or {}).get("value")
            if value:
                return value.strip()
    return ""


def entries(accessions: list[str], batch: int = 80, max_age_s: float = 30 * 86400) -> dict[str, dict | None]:
    """{accession: {"name": UniProt's protein name, "seq": sequence of the current entry}}.

    None means UniProt has no such entry any more (deleted or merged away). That is not
    the same as an entry without membrane segments, which membrane_features cannot tell
    apart: it answers {"tm": [], "named": []} for both."""
    path = CACHE_DIR / "uniprot_entry.json"
    cache: dict = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    todo = sorted({a for a in accessions if a and re.fullmatch(r"[A-Z0-9]{6,10}", a)} - cache.keys())
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        url = ACCESSIONS + "?" + urllib.parse.urlencode({"accessions": ",".join(chunk), "format": "json", "size": 500,
                                                         "fields": "accession,protein_name,sequence"})
        results = json.loads(_get(url)).get("results") or []
        for a in chunk:
            cache[a] = None
        for x in results:
            rec = {"name": _protein_name(x.get("proteinDescription")), "seq": (x.get("sequence") or {}).get("value") or ""}
            if not rec["seq"]:  # an inactive entry
                continue
            for a in [x.get("primaryAccession")] + (x.get("secondaryAccessions") or []):
                if a in chunk and (cache.get(a) is None or a == x.get("primaryAccession")):
                    cache[a] = rec
        if (i // batch) % 20 == 0 or i + batch >= len(todo):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    return {a: cache.get(a) for a in accessions}
