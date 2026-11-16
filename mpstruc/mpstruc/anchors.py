"""Links to listings on the official mpstruc site.

The mpstruc home page gives every parent and daughter listing a bookmark
anchor, id_<PDB>_<subgroup number>_<group number> (the page's own bookmark
icons link to https://blanco.biomol.uci.edu/mpstruc/#id_1XKW_14_2), and every
subfamily an anchor id_<GROUP>:<Subfamily>. Related structures have no anchor
of their own: they are a line in their owner's listing.

    python -m mpstruc.anchors [page.html]   # write data/mpstruc_links.json

With no argument the page is downloaded. The weekly job runs this; if the
site cannot be reached the existing file is kept.
"""
from __future__ import annotations

import json
import re
import sys
import urllib.request
from pathlib import Path

HOME = "https://blanco.biomol.uci.edu/mpstruc/"
OUT = Path(__file__).resolve().parents[2] / "data" / "mpstruc_links.json"
_ENTRY = re.compile(r'id="id_([0-9][A-Z0-9]{3})_(\d+)_(\d+)"')
_SECTION = re.compile(r'id="(id_[A-Z_]+:[^"]+)"')


def parse(html: str) -> dict:
    """{"base": HOME, "entries": {PDB: "76_2"}, "sections": [anchor ids]} from the page."""
    entries = {pdb: f"{sg}_{g}" for pdb, sg, g in _ENTRY.findall(html)}
    return {"base": HOME, "entries": entries, "sections": sorted(set(_SECTION.findall(html)))}


def main() -> None:
    if len(sys.argv) > 1:
        html = Path(sys.argv[1]).read_text(encoding="utf-8", errors="replace")
    else:
        req = urllib.request.Request(HOME, headers={"User-Agent": "mpstruc-curation-tool"})
        with urllib.request.urlopen(req, timeout=300) as r:
            html = r.read().decode("utf-8", "replace")
    links = parse(html)
    if len(links["entries"]) < 1000:  # a broken or partial page: keep what we have
        raise SystemExit(f"only {len(links['entries'])} listing anchors found; keeping the existing file")
    OUT.write_text(json.dumps(links, separators=(",", ":")))
    print(f"wrote {OUT} ({len(links['entries'])} listings, {len(links['sections'])} sections)")


if __name__ == "__main__":
    main()
