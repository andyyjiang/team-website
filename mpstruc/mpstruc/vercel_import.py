"""Keep a local snapshot of the other group's Vercel demo (mpstruc-viewer-demo).

    python -m mpstruc.vercel_import [export.csv] [--offline]
    python -m mpstruc.vercel_import page <groups.jsonl | page.txt> [--offline]

The first form imports a CSV export of their database (default: the newest
file in vercel/exports/). The second merges what vercel/read_page.js read
from the rendered demo page (JSON lines), or a plain-text copy of the page.
Both update ../data/vercel_snapshot.json:

    {"source": "export:<file>" | "page", "readAt": ISO, "export": <file>|null,
     "pageReadAt": ISO|null,
     "entries": {PDB: {"status", "family", "name", "group", "released", "tm",
                       "memo", "seenAt"}}}

An export replaces the fields it has (status, family, name, group, memo) and
release dates are filled from RCSB; a page read adds entries and updates
statuses on top of it. seenAt is the first date we saw a PDB listed.

The Vercel curators' memos are the other group's working notes and this
repository is public, so they are never written to data/vercel_snapshot.json
(memo is null there); they are kept in cache/vercel_memos.json (not committed)
and put back when the snapshot is loaded. The CSV exports themselves stay in
vercel/exports/, which git ignores.

A page read with no groups changes nothing and exits non-zero.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import sys
from pathlib import Path

from . import rcsb

ROOT = Path(__file__).resolve().parent.parent
EXPORTS = ROOT / "vercel" / "exports"
SNAPSHOT = ROOT.parent / "data" / "vercel_snapshot.json"
MEMOS = rcsb.CACHE_DIR / "vercel_memos.json"

STATUS_CODES = {
    "already in": "already_in",
    "no": "no",
    "yes": "yes",
    "maybe": "maybe",
    "pubmed ready": "pubmed_ready",
    "ready for yes": "ready_for_yes",
    "": "pending",
    "(pending)": "pending",
    "pending": "pending",
}
PAGE_VERDICTS = {"Yes": "yes", "No": "no", "Maybe": "maybe", "(Pending)": "pending"}
# Export statuses that say more than the page's verdict of the same direction:
# the page shows (Pending) or Yes for them, so a page read must not overwrite them.
_REFINES = {"pending": {"already_in", "pubmed_ready", "ready_for_yes", "yes", "no", "maybe"},
            "yes": {"already_in", "ready_for_yes"}}
_PDB = re.compile(r"^[0-9][A-Z0-9]{3}$|^PDB_[0-9A-Z]{8}$|^[0-9]{2}[A-Z0-9]{2}$")


def status_code(raw: str | None) -> str:
    s = (raw or "").strip().lower()
    return STATUS_CODES.get(s) or re.sub(r"[^a-z0-9]+", "_", s).strip("_")


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def _blank(v: str | None) -> str | None:
    v = (v or "").strip()
    return v or None


# --------------------------------------------------------------------------- RCSB

def release_info(pdb_ids, batch: int = 300, offline: bool = False) -> dict[str, dict | None]:
    """{PDB: {"released": "YYYY-MM-DD", "title": str} | None} from RCSB GraphQL.

    Cached in cache/release_dates.json (release dates never change, so found
    entries are kept for good; unknown ones are asked again next time). The
    entry_summaries cache seeds it so already-known entries cost no request.
    """
    path = rcsb.CACHE_DIR / "release_dates.json"
    cache: dict[str, dict | None] = json.loads(path.read_text()) if path.exists() else {}
    ids = sorted({p.upper() for p in pdb_ids})
    missing = [p for p in ids if not cache.get(p)]
    if missing:
        seed = rcsb.CACHE_DIR / "entry_summaries.json"
        if seed.exists():
            summ = json.loads(seed.read_text())
            for p in missing:
                e = summ.get(p)
                rel = ((e or {}).get("rcsb_accession_info") or {}).get("initial_release_date")
                if rel:
                    cache[p] = {"released": rel[:10], "title": ((e.get("struct") or {}).get("title"))}
        missing = [p for p in ids if not cache.get(p)]
    if missing and not offline:
        q = ("query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id struct { title } "
             "rcsb_accession_info { initial_release_date } } }")
        for i in range(0, len(missing), batch):
            chunk = missing[i:i + batch]
            data = rcsb._post_json(rcsb.GRAPHQL, {"query": q, "variables": {"ids": chunk}})
            for e in (data.get("data") or {}).get("entries") or []:
                if not e:
                    continue
                rel = (e.get("rcsb_accession_info") or {}).get("initial_release_date")
                cache[e["rcsb_id"]] = {"released": rel[:10] if rel else None,
                                       "title": (e.get("struct") or {}).get("title")}
            for p in chunk:
                cache.setdefault(p, None)
            print(f"  RCSB release dates {min(i + batch, len(missing))}/{len(missing)}", file=sys.stderr)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cache, separators=(",", ":")))
    return {p: cache.get(p) for p in ids}


# --------------------------------------------------------------------------- snapshot

def _memo_path(path: Path) -> Path:
    """Where the memos of the snapshot at path are kept (outside the published data)."""
    return MEMOS if path == SNAPSHOT else path.with_name(path.stem + ".memos.json")


def load_snapshot(path: Path = SNAPSHOT, memos: Path | None = None) -> dict:
    if path.exists():
        snap = json.loads(path.read_text())
        snap.setdefault("entries", {})
        mp = memos or _memo_path(path)
        kept = json.loads(mp.read_text()) if mp.exists() else {}
        for pdb, memo in kept.items():
            if pdb in snap["entries"] and not snap["entries"][pdb].get("memo"):
                snap["entries"][pdb]["memo"] = memo
        return snap
    return {"source": None, "readAt": None, "export": None, "pageReadAt": None, "entries": {}}


def save_snapshot(snap: dict, path: Path = SNAPSHOT, memos: Path | None = None) -> None:
    """Write the snapshot without the Vercel curators' memos (kept in memos, not published)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    snap["entries"] = dict(sorted(snap["entries"].items()))
    kept = {p: e["memo"] for p, e in snap["entries"].items() if e.get("memo")}
    mp = memos or _memo_path(path)
    mp.parent.mkdir(parents=True, exist_ok=True)
    mp.write_text(json.dumps(kept, ensure_ascii=False, separators=(",", ":")))
    public = {**snap, "entries": {p: {**e, "memo": None} for p, e in snap["entries"].items()}}
    path.write_text(json.dumps(public, ensure_ascii=False, separators=(",", ":")) + "\n")


def _entry(snap: dict, pdb: str) -> dict:
    return snap["entries"].setdefault(pdb, {"status": "pending", "family": None, "name": None, "group": None,
                                            "released": None, "tm": None, "memo": None, "seenAt": None})


def _seen(e: dict, day: str) -> None:
    if not e.get("seenAt") or day < e["seenAt"]:
        e["seenAt"] = day


def counts(snap: dict) -> dict:
    by_status: dict[str, int] = {}
    for e in snap["entries"].values():
        by_status[e["status"]] = by_status.get(e["status"], 0) + 1
    return {"entries": len(snap["entries"]), "byStatus": dict(sorted(by_status.items())),
            "withReleaseDate": sum(1 for e in snap["entries"].values() if e.get("released"))}


# --------------------------------------------------------------------------- CSV export

def newest_export(folder: Path = EXPORTS) -> Path:
    files = sorted(folder.glob("*.csv"), key=lambda p: (_file_date(p) or "", p.stat().st_mtime))
    if not files:
        raise SystemExit(f"no CSV export in {folder}")
    return files[-1]


def _file_date(p: Path) -> str | None:
    m = re.search(r"(\d{4}-\d{2}-\d{2})", p.name)
    return m.group(1) if m else None


def _pmid(row: dict) -> str | None:
    m = re.search(r"PMID:\s*(\d+)", row.get("Bibliography") or "")
    if not m:
        m = re.search(r"[?&]term=(\d+)", row.get("PubMed Link") or "")
    return m.group(1) if m else None


def read_export(path: Path) -> dict[str, dict]:
    """{PDB: {status, family, name, group, memo}} from a Vercel CSV export."""
    out: dict[str, dict] = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            pdb = (row.get("Pdb Code") or "").strip().upper()
            if not pdb:
                continue
            out[pdb] = {"status": status_code(row.get("Status")), "family": _blank(row.get("Family")),
                        "name": _blank(row.get("Name")), "group": _pmid(row), "memo": _blank(row.get("Memo"))}
    return out


def import_export(path: Path, snap: dict, offline: bool = False) -> dict:
    rows = read_export(path)
    day = _file_date(path) or dt.date.today().isoformat()
    info = release_info(rows, offline=offline)
    for pdb, r in rows.items():
        e = _entry(snap, pdb)
        e.update(r)
        rel = (info.get(pdb) or {}).get("released")
        if rel:
            e["released"] = rel
        _seen(e, day)
    snap["source"] = f"export:{path.name}"
    snap["export"] = path.name
    snap["readAt"] = _now()
    return {"rows": len(rows), "noReleaseDate": sorted(p for p in rows if not snap["entries"][p].get("released"))}


# --------------------------------------------------------------------------- rendered page

_FIELD = re.compile(r"^(Title|Release Date\(s\)|Number of TM segment|Resolution|UniProt ID\(s\)|Classification\(s\)"
                    r"|Membrane Protein Verdict):\s*(.*)$")
_PRED = re.compile(r"^(.+? model on [^:]+):\s*(.*?)\s+[—-]\s+([\d.]+)%$")
_SPLIT = re.compile(r",\s*")


def parse_group(text: str) -> dict | None:
    """One group card's visible text -> dict (same shape read_page.js emits)."""
    lines = [ln.strip() for ln in text.splitlines()]
    lines = [ln for ln in lines if ln]
    # Lines above the label (a "#1" badge, say) are not part of the group.
    at = next((i for i, ln in enumerate(lines) if ln.startswith("PubMed Group:")), None)
    if at is None:
        return None
    lines = lines[at:]
    g = {"group": _blank(lines[0].split(":", 1)[1]), "journal": None, "paperTitle": None, "doi": None,
         "title": None, "pdbIds": [], "releaseDates": [], "tm": None, "resolution": None, "uniprot": [],
         "classifications": [], "predictions": [], "verdict": None}
    i, n = 1, len(lines)
    head: list[str] = []
    while i < n and not lines[i].startswith("Title:") and lines[i] != "PDB IDs:":
        head.append(lines[i])
        i += 1
    if head and head[0] == "Journal:":
        head = head[1:]
        if head:
            g["journal"] = None if head[0] in ("N/A", "") else head[0]
            head = head[1:]
    urls = [h for h in head if h.startswith("http")]
    g["doi"] = urls[0] if urls else None
    g["paperTitle"] = " ".join(h for h in head if not h.startswith("http")) or None
    while i < n:
        ln = lines[i]
        if ln == "PDB IDs:":
            i += 1
            while i < n and _PDB.match(lines[i].upper()):
                g["pdbIds"].append(lines[i].upper())
                i += 1
            continue
        m = _FIELD.match(ln)
        if m:
            key, val = m.group(1), m.group(2).strip()
            if key == "Title":
                g["title"] = val or None
            elif key == "Release Date(s)":
                g["releaseDates"] = re.findall(r"\d{4}-\d{2}-\d{2}", val)
            elif key == "Number of TM segment":
                g["tm"] = val or None
            elif key == "Resolution":
                g["resolution"] = val or None
            elif key == "UniProt ID(s)":
                g["uniprot"] = [u for u in _SPLIT.split(val) if u and u != "N/A"]
            elif key == "Classification(s)":
                g["classifications"] = [c for c in _SPLIT.split(val) if c]
            else:
                g["verdict"] = val or None
        else:
            p = _PRED.match(ln)
            if p:
                g["predictions"].append({"model": p.group(1), "label": p.group(2), "pct": float(p.group(3))})
        i += 1
    return g


def parse_page_text(text: str) -> list[dict]:
    """Every group in a plain-text copy of the rendered page (document.body.innerText)."""
    parts = re.split(r"(?m)^(?=\s*PubMed Group:)", text)
    return [g for g in (parse_group(p) for p in parts) if g]


def load_page_file(path: Path) -> tuple[list[dict], dict]:
    """Groups and the reader's meta line from read_page.js output, or from a page text copy."""
    text = Path(path).read_text(encoding="utf-8")
    if text.lstrip().startswith("{"):
        groups, meta = [], {}
        for ln in text.splitlines():
            if not ln.strip():
                continue
            obj = json.loads(ln)
            if "_meta" in obj:
                meta = obj["_meta"]
            else:
                groups.append(obj)
        return groups, meta
    return parse_page_text(text), {}


def merge_page(groups: list[dict], snap: dict, seen_day: str | None = None, offline: bool = False) -> dict:
    """Add/update snapshot entries from page groups; returns counts.

    No groups (the page could not be read) is an error: the snapshot is left as it
    was, so it never claims a page read that did not happen."""
    if not any(g.get("pdbIds") for g in groups):
        raise ValueError("the page read holds no group with PDB IDs; the snapshot is unchanged")
    day = seen_day or dt.date.today().isoformat()
    ids = {p for g in groups for p in g.get("pdbIds") or []}
    info = release_info(ids, offline=offline)
    added = changed = 0
    for g in groups:
        code = PAGE_VERDICTS.get((g.get("verdict") or "").strip(), status_code(g.get("verdict")) or "pending")
        pdbs = g.get("pdbIds") or []
        dates = g.get("releaseDates") or []
        classes = sorted({c for c in g.get("classifications") or [] if c != "N/A"})
        group = g.get("group") if re.fullmatch(r"\d+", g.get("group") or "") else None
        for pdb in pdbs:
            new = pdb not in snap["entries"]
            e = _entry(snap, pdb)
            added += new
            old = e["status"]
            if new or old not in _REFINES.get(code, set()):
                e["status"] = code
            changed += (not new and e["status"] != old)
            if group:
                e["group"] = group
            if g.get("tm"):
                e["tm"] = g["tm"]
            rel = (info.get(pdb) or {}).get("released") or (dates[0] if len(set(dates)) == 1 else None)
            if rel and not e.get("released"):
                e["released"] = rel
            if not e.get("name"):
                e["name"] = (info.get(pdb) or {}).get("title") or (g.get("title") if len(pdbs) == 1 else None)
            if not e.get("family") and len(classes) == 1:
                e["family"] = classes[0]
            _seen(e, day)
    snap["source"] = "page"
    snap["pageReadAt"] = _now()
    snap["readAt"] = snap["pageReadAt"]
    return {"groups": len(groups), "pdbs": len(ids), "added": added, "statusChanged": changed}


# --------------------------------------------------------------------------- CLI

def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("args", nargs="*", help="[export.csv] or: page <file>")
    ap.add_argument("--offline", action="store_true", help="do not ask RCSB for missing release dates")
    ap.add_argument("--snapshot", type=Path, default=SNAPSHOT)
    a = ap.parse_args(argv)
    snap = load_snapshot(a.snapshot)
    if a.args and a.args[0] == "page":
        if len(a.args) != 2:
            ap.error("usage: page <file>")
        groups, meta = load_page_file(Path(a.args[1]))
        seen = (meta.get("readAt") or "")[:10] or None
        try:
            res = merge_page(groups, snap, seen_day=seen, offline=a.offline)
        except ValueError as e:
            raise SystemExit(f"vercel_import: {e} ({meta.get('stopReason') or 'no reader notes'})")
        if meta:
            res["pageTotalGroups"] = meta.get("total")
    else:
        path = Path(a.args[0]) if a.args else newest_export()
        res = import_export(path, snap, offline=a.offline)
        res["noReleaseDate"] = len(res["noReleaseDate"])
    save_snapshot(snap, a.snapshot)
    res["snapshot"] = counts(snap)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
