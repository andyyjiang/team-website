"""Automatic "Not published yet" rejections by the weekly run.

mpstruc lists a structure only once its paper is out. The weekly run
(automation.main) finds which suggestions are still waiting for their paper
(awaiting.py). Instead of leaving them in the curators' suggestion list,
this module records them in data/curation.json as rejected by the
automation:

    entries[PDB] = {"status": "rejected", "reasonCode": "unpublished",
                    "reason": "Not published yet (PDB lists it as 'To be published')",
                    "by": "mpstruc-automation", "byName": "Automation (Claude)",
                    "at": ISO time, "auto": true, ...snapshot of the suggestion}

and brings them back when the paper is published: a rejection with
reasonCode "unpublished" is removed from data/curation.json once awaiting.py
sees the paper appear, so the structure returns to "To decide" with its
"Paper now published" badge. A curator's own "Not published yet" rejection
is lifted only if the run saw the structure still waiting after the
curator's decision and the paper appeared later (so a curator who rejects a
structure whose PDB entry already names a journal, e.g. one without page
numbers yet, is not overruled every week). Approved entries, deleted ones
(tombstones) and rejections for any other reason are never touched, and a
structure a curator moved back to "To decide" (unarchive, restore from
Deleted, or undo) is not rejected again.

The site may commit curator decisions to data/curation.json while the job
runs, so the file is changed by read-modify-write: only the automation's own
entries are added or removed and log lines appended; everything else is
written back as read (the site's own format: JSON with indent 1). If the
job's push conflicts on that file, `publish` takes the newer file from the
repository and re-applies the change (`python -m mpstruc.autoreject apply`
does the same by hand from data/candidates.json and data/awaiting.json).

    python -m mpstruc.autoreject apply              # update data/curation.json
    python -m mpstruc.autoreject push --branch main # commit the data files and push

The functions above the CLI section are pure (no files, no network).
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Callable

from . import awaiting

AUTO_BY = "mpstruc-automation"
AUTO_NAME = "Automation (Claude)"
UNPUBLISHED = "unpublished"
LABEL = "Not published yet"
CURATION = "data/curation.json"
# Data files the weekly job writes and commits (paths from the repository root).
DATA_FILES = ["data/candidates.json", "data/mpstruc.json", "data/import_log.json",
              "data/awaiting.json", CURATION]


def now_iso() -> str:
    """The time the way the site writes it (2026-11-06T00:11:16.810Z)."""
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def reason_text(why: str | None) -> str:
    if why == awaiting.PREPRINT_ONLY:
        return f"{LABEL} (preprint only)"
    return f"{LABEL} (PDB lists it as 'To be published')"


def _why_phrase(why: str | None) -> str:
    if why == awaiting.PREPRINT_ONLY:
        return "not published yet (only a preprint so far)"
    return "not published yet (PDB lists it as 'To be published')"


def is_unpublished_rejection(e) -> bool:
    return isinstance(e, dict) and e.get("status") == "rejected" and e.get("reasonCode") == UNPUBLISHED


def unpublished_rejections(curation: dict | None) -> dict[str, dict]:
    """Rejections waiting for a paper (automatic or by a curator), by PDB ID."""
    ents = (curation or {}).get("entries") or {}
    return {k.upper(): v for k, v in ents.items() if is_unpublished_rejection(v)} if isinstance(ents, dict) else {}


def decided_ids(curation: dict | None) -> set[str]:
    """Structures the run must not suggest or score: everything in
    data/curation.json except rejections waiting for a paper, which are
    still checked every week so they can come back when it is out."""
    ents = (curation or {}).get("entries") or {}
    if not isinstance(ents, dict):
        return set()
    return {k.upper() for k, v in ents.items() if not is_unpublished_rejection(v)}


def track_unpublished(rec: dict, curation: dict | None, today: str) -> tuple[dict, list[str]]:
    """Add rejections waiting for a paper to the awaiting record (awaiting.py)
    if it does not track them yet (e.g. a curator rejected one as "Not
    published yet"), so later runs re-check them. Structures released (or
    rejected) longer ago than awaiting.GIVE_UP_DAYS are not added again.

    Added entries are marked "unconfirmed": the run has not yet seen these
    structures waiting for a paper. If this run finds the PDB citation
    already counts as published, awaiting.update just drops them again
    (no "newly published" date), so the curator's rejection stays."""
    cutoff = (dt.date.fromisoformat(today) - dt.timedelta(days=awaiting.GIVE_UP_DAYS)).isoformat()
    wait = dict(rec.get("awaiting") or {})
    pub = rec.get("published") or {}
    added = []
    for pdb, e in sorted(unpublished_rejections(curation).items()):
        if pdb in wait or pdb in pub:
            continue
        start = str(e.get("released") or e.get("at") or "")[:10]
        if start and start < cutoff:
            continue
        why = awaiting.PREPRINT_ONLY if "preprint" in str(e.get("reason") or "").lower() else awaiting.TO_BE_PUBLISHED
        wait[pdb] = {"since": str(e.get("at") or today)[:10], "released": e.get("released"), "reason": why,
                     "unconfirmed": True}
        added.append(pdb)
    return {**rec, "awaiting": wait}, added


def snapshot(c: dict) -> dict:
    """The suggestion as saved with a decision (same fields as the site's snapshot())."""
    return {"title": c.get("title"), "released": c.get("released"), "organism": (c.get("organisms") or [""])[0] or "",
            "res": c.get("res"), "method": c.get("method"), "pmid": c.get("pmid"), "doi": c.get("doi"),
            "mp": c.get("mp"), "suggested": c.get("placement"), "suggestedSubgroup": c.get("subgroup"),
            "suggestedOwner": c.get("ownerPdb") or "", "membrane": c.get("membrane")}


# Curator actions that put a structure back in "To decide" (when it has no
# entry in data/curation.json afterwards).
_BACK = ("unarchive", "restore")


def curator_unarchived(log: list, pdb: str) -> bool:
    """True if a curator moved this structure back to "To decide" after the
    automation last rejected it (or without it ever being rejected
    automatically): the automation then leaves it alone.

    Called only for structures with no entry now. "unarchive" and "restore"
    (from Deleted) always count; "undo" counts when it came after an
    automatic rejection (e.g. a curator undid the deletion of one)."""
    back = None
    for line in reversed(log if isinstance(log, list) else []):
        if not isinstance(line, dict) or str(line.get("pdb") or "").upper() != pdb:
            continue
        if line.get("by") == AUTO_BY:
            if line.get("action") == "auto-reject":
                return back is not None
            continue
        if back is None and line.get("action") in (*_BACK, "undo"):
            back = line.get("action")
    return back in _BACK


def _lift(e: dict, on: str) -> bool:
    """Whether a "Not published yet" rejection should be lifted now that the
    run saw the paper on date `on`. The automation's own rejections always
    are; a curator's only if the paper appeared after their decision (so a
    curator who rejects it again after it came back is not overruled)."""
    if e.get("auto") and e.get("by") == AUTO_BY:
        return True
    return on > str(e.get("at") or "")[:10]


def published_on(pdb: str, cand: dict | None, rec: dict) -> str | None:
    """The date the run first saw this structure's paper, or None if it is
    still waiting (or nothing is known)."""
    if cand is not None:
        if cand.get("awaiting"):
            return None
        if cand.get("publishedOn"):
            return cand["publishedOn"]
    if pdb in (rec.get("awaiting") or {}):
        return None
    on = ((rec.get("published") or {}).get(pdb) or {}).get("on")
    return on or None


def update_curation(raw: dict, cands: list[dict], rec: dict, at: str) -> dict:
    """Add and remove the automation's "Not published yet" rejections in a
    freshly read data/curation.json (raw, changed in place).

    cands: this run's candidates (with "awaiting" / "publishedOn" from
    awaiting.update); rec: the awaiting record after this run.
    Returns {"rejected": [...], "unrejected": [...]} (PDB IDs).
    Running it again on its own output changes nothing.
    """
    if not isinstance(raw.get("entries"), dict):
        raw["entries"] = {}
    if not isinstance(raw.get("log"), list):
        raw["log"] = []
    ents, log = raw["entries"], raw["log"]
    keys = {k.upper(): k for k in ents}
    by_pdb = {c["pdb"].upper(): c for c in cands if c.get("pdb")}
    rejected, unrejected = [], []
    # 1. Suggestions still waiting for their paper, not decided yet.
    for pdb, c in sorted(by_pdb.items()):
        wait = c.get("awaiting")
        if not wait or pdb in keys or curator_unarchived(log, pdb):
            continue
        why = wait.get("reason")
        ents[pdb] = {"status": "rejected", "reasonCode": UNPUBLISHED, "reason": reason_text(why),
                     "by": AUTO_BY, "byName": AUTO_NAME, "at": at, "auto": True, **snapshot(c)}
        log.append({"at": at, "by": AUTO_BY, "byName": AUTO_NAME, "action": "auto-reject", "pdb": pdb,
                    "detail": f"Rejected {pdb} automatically: {_why_phrase(why)}. "
                              "It will come back to To decide when the paper is out."})
        rejected.append(pdb)
    # 2. Rejections waiting for a paper whose paper is now out.
    for pdb, key in sorted(keys.items()):
        if not is_unpublished_rejection(ents.get(key)):
            continue
        on = published_on(pdb, by_pdb.get(pdb), rec)
        if not on or not _lift(ents[key], on):
            continue
        del ents[key]
        log.append({"at": at, "by": AUTO_BY, "byName": AUTO_NAME, "action": "auto-unreject", "pdb": pdb,
                    "detail": f"Paper for {pdb} is now published ({on}): back in To decide"})
        unrejected.append(pdb)
    return {"rejected": rejected, "unrejected": unrejected}


def mark_candidates(cands: list[dict], raw: dict) -> None:
    """Set c["autoRejected"] on candidates the automation has rejected as not
    published yet (the site shows them in the Rejected tab)."""
    ents = raw.get("entries") or {}
    up = {k.upper(): v for k, v in ents.items()} if isinstance(ents, dict) else {}
    for c in cands:
        e = up.get(str(c.get("pdb") or "").upper())
        if is_unpublished_rejection(e) and e.get("auto"):
            c["autoRejected"] = True
        else:
            c.pop("autoRejected", None)


_SURROGATE = re.compile("[\ud800-\udfff]")


def dumps(raw: dict) -> str:
    """data/curation.json as the site writes it: JSON.stringify(raw, null, 1) + "\\n".

    Like JSON.stringify, a lone surrogate (e.g. half an emoji cut off by a
    text field's length limit) is written as a \\uXXXX escape: it cannot be
    written as UTF-8. (Python never holds a valid pair as two surrogates.)"""
    text = json.dumps(raw, indent=1, ensure_ascii=False)
    return _SURROGATE.sub(lambda m: "\\u%04x" % ord(m.group()), text) + "\n"


def update_text(text: str | None, cands: list[dict], rec: dict, at: str) -> tuple[str | None, dict, dict]:
    """update_curation on the file's text. Returns (new text, or the input
    unchanged when nothing changed; the parsed result; counts).
    Raises ValueError if the file is not the expected JSON object."""
    if text is None or not text.strip():
        raw, text = {"version": 1, "entries": {}, "log": []}, None
    else:
        raw = json.loads(text)
        if not isinstance(raw, dict) or ("entries" in raw and not isinstance(raw["entries"], dict)) \
                or ("log" in raw and not isinstance(raw["log"], list)):
            raise ValueError("data/curation.json does not have the expected shape")
    counts = update_curation(raw, cands, rec, at)
    if not counts["rejected"] and not counts["unrejected"]:
        return text, raw, counts
    return dumps(raw), raw, counts


def set_stats(stats: dict, counts: dict) -> None:
    """This run's counts in data/candidates.json's "stats"."""
    stats["rejected automatically: not published yet"] = len(counts["rejected"])
    stats["back after publication"] = len(counts["unrejected"])


def commit_message(day: str, n_rejected: int, n_back: int) -> str:
    return (f"Weekly mpstruc update ({day}): {n_rejected} auto-rejected as not published yet, "
            f"{n_back} back after publication")


# ------------------------------------------------------------- git publishing

def _git(repo: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True)


def publish(repo: Path, branch: str, files: list[str], message: Callable[[], str],
            reapply: Callable[[], None], attempts: int = 3, remote: str = "origin") -> str:
    """Commit the data files and push them to branch.

    Curators may save decisions to data/curation.json while the job runs.
    The commit is rebased onto them; if that conflicts on data/curation.json,
    the rebase is aborted, the commit is rebuilt on the newest branch with
    that file taken from the repository and reapply() run again (it re-adds
    only the automation's own changes), up to `attempts` times.
    Returns "pushed", "no changes", or raises RuntimeError.
    """
    def ok(*args):
        r = _git(repo, *args)
        if r.returncode:
            raise RuntimeError(f"git {' '.join(args)} failed: {r.stderr.strip() or r.stdout.strip()}")
        return r.stdout

    def stage():
        ok("add", "--", *[f for f in files if (Path(repo) / f).exists()])
        return _git(repo, "diff", "--cached", "--quiet").returncode != 0

    if not stage():
        return "no changes"
    ok("commit", "-q", "-m", message())
    for _ in range(attempts):
        pulled = _git(repo, "pull", "-q", "--rebase", remote, branch)
        if pulled.returncode:
            conflicted = set(_git(repo, "diff", "--name-only", "--diff-filter=U").stdout.split())
            _git(repo, "rebase", "--abort")
            if conflicted - {CURATION}:
                raise RuntimeError(f"git pull --rebase failed: {', '.join(sorted(conflicted)) or pulled.stderr.strip()}")
            # Rebuild the commit on the newest branch: its data/curation.json
            # plus the automation's own changes, and this run's other files.
            ok("fetch", "-q", remote, branch)
            ok("reset", "-q", "--mixed", "FETCH_HEAD")
            ok("checkout", "FETCH_HEAD", "--", CURATION)
            reapply()
            if not stage():
                return "no changes"
            ok("commit", "-q", "-m", message())
        if _git(repo, "push", "-q", remote, f"HEAD:{branch}").returncode == 0:
            return "pushed"
    raise RuntimeError(f"could not push to {branch} after {attempts} attempts")


# ---------------------------------------------------------------------- CLI

def _data_dir() -> Path:
    from .report import ROOT
    return ROOT.parent / "data"


def apply(data: Path) -> dict:
    """Re-apply this run's changes to data/curation.json from
    data/candidates.json and data/awaiting.json. Also updates the candidates'
    autoRejected marks and this run's counts in data/import_log.json."""
    cand_path, cur_path = data / "candidates.json", data / "curation.json"
    out = json.loads(cand_path.read_text(encoding="utf-8"))
    cands = out.get("candidates") or []
    rec = awaiting.load(data / "awaiting.json")
    text = cur_path.read_text(encoding="utf-8") if cur_path.exists() else None
    new, raw, counts = update_text(text, cands, rec, now_iso())
    if new != text:
        cur_path.write_text(new, encoding="utf-8")
    mark_candidates(cands, raw)
    if isinstance(out.get("stats"), dict):
        set_stats(out["stats"], counts)
    cand_path.write_text(json.dumps(out, separators=(",", ":")), encoding="utf-8")
    log_path = data / "import_log.json"
    if log_path.exists():
        runs = json.loads(log_path.read_text()).get("runs", [])
        if runs and runs[-1].get("at") == out.get("generated"):
            runs[-1]["autoRejected"] = len(counts["rejected"])
            runs[-1]["autoUnrejected"] = len(counts["unrejected"])
            log_path.write_text(json.dumps({"runs": runs}, indent=1))
    return counts


def _last_counts(data: Path) -> tuple[int, int]:
    try:
        run = (json.loads((data / "import_log.json").read_text()).get("runs") or [{}])[-1]
    except (OSError, ValueError):
        run = {}
    return int(run.get("autoRejected") or 0), int(run.get("autoUnrejected") or 0)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m mpstruc.autoreject")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("apply", help="update data/curation.json from data/candidates.json and data/awaiting.json")
    p = sub.add_parser("push", help="commit the weekly data files and push them")
    p.add_argument("--branch", required=True)
    p.add_argument("--attempts", type=int, default=3)
    args = ap.parse_args(argv)
    data = _data_dir()
    if args.cmd == "apply":
        counts = apply(data)
        print(f"curation: {len(counts['rejected'])} rejected as not published yet, "
              f"{len(counts['unrejected'])} back after publication")
        return
    day = dt.datetime.now(dt.timezone.utc).date().isoformat()
    result = publish(data.parent, args.branch, DATA_FILES, lambda: commit_message(day, *_last_counts(data)),
                     lambda: apply(data), attempts=args.attempts)
    print(result)


if __name__ == "__main__":
    main(sys.argv[1:])
