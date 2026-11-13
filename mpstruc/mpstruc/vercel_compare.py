"""Compare our membrane-protein evidence with the Vercel demo's list.

    python -m mpstruc.vercel_compare                    # the latest PDB release week
    python -m mpstruc.vercel_compare --week 2026-11-04  # a given release Wednesday
    python -m mpstruc.vercel_compare --whole            # Vercel's whole snapshot (long first run)
    python -m mpstruc.vercel_compare --week 2026-11-11 --whole   # both

Reads data/vercel_snapshot.json (written by the Vercel importer; when it is
missing, the newest export in mpstruc/vercel/exports/ is read directly),
data/candidates.json, data/curation.json, the official mpstruc snapshot
(report.source_xml), PDBTM and MemProtMD (through the PDB's search service),
OPM (membranedb.opm) and our membrane model (the same steps as
membrane.score_structures, with one context for all batches). Writes
data/vercel_compare.json. Past weeks are kept; the week being compared is
replaced, and the whole-database section is kept unless --whole runs again.

The curators' rule: a structure is a membrane protein only when a significant
transmembrane segment(s) or membrane-spanning beta-barrel is present in it.
Each structure's verdict comes from the strongest evidence there is, in this
order:

1. our curators' decision (approved => membrane; rejected as not a membrane
   protein => not a membrane protein; other rejections, such as "not published
   yet", say nothing about the membrane and are shown only as evidence);
2. listed in the official mpstruc => membrane, unless the structure at that
   code has no resolved membrane-crossing segment or barrel (a soluble fragment
   mpstruc lists, or a typo in mpstruc's PDB code): then mpstruc gives no vote;
3. PDBTM, MemProtMD, or OPM placing it across the membrane => membrane (OPM
   only when it finds at least one membrane-crossing segment in the structure;
   MemProtMD alone, with nothing membrane-crossing resolved and our model under
   10%, makes it uncertain);
4. our model: 50% or more => membrane, under 10% => not a membrane protein
   (with the resolved-segment facts as reasons), else uncertain.

"Missed" = verdict membrane, released in the compared period, and nowhere in
Vercel's snapshot. "notMembrane" = in Vercel's snapshot with verdict "not a
membrane protein" (with Vercel's own status, so readers see whether Vercel
already said No). "agree" (per week) = structures of the week Vercel lists that
our evidence also calls membrane proteins.

Whole-database mode checks every structure in Vercel's snapshot (scoring those
no stronger source decides; model scores are checkpointed to
cache/vercel_compare_scores.json, so an interrupted run resumes), and looks for
missed structures among those released within the span Vercel covers using the
strong sources only (official mpstruc, PDBTM, MemProtMD, OPM transmembrane),
plus our model for the last 26 weeks of that span, so the whole PDB is not
scored.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from .rcsb import CACHE_DIR
from .report import ROOT

SITE_DATA = ROOT.parent / "data"
SNAPSHOT = SITE_DATA / "vercel_snapshot.json"
OUTPUT = SITE_DATA / "vercel_compare.json"
EXPORTS = ROOT / "vercel" / "exports"
SCORE_CACHE = CACHE_DIR / "vercel_compare_scores.json"
ANNOT_CACHE = CACHE_DIR / "vercel_compare_annotated.json"

MEMBRANE, NOT_MEMBRANE, UNCERTAIN = "membrane", "not_membrane", "uncertain"
P_MEMBRANE, P_NOT = 0.5, 0.1
RECENT_WEEKS = 26
STATUS_FROM_EXPORT = {"already in": "already_in", "no": "no", "yes": "yes", "maybe": "maybe",
                      "pubmed ready": "pubmed_ready", "ready for yes": "ready_for_yes", "": "pending"}
OPM_TEXT = {
    "transmembrane": "OPM places it across the membrane (transmembrane)",
    "transmembrane-family, no segment in structure":
        "OPM files it with transmembrane proteins but finds no membrane-crossing segment in this structure, "
        "so OPM gives no vote",
    "monotopic/peripheral": "OPM places it only on the membrane surface (monotopic/peripheral), "
                            "which on its own does not meet the curators' rule",
    "peptide": "OPM lists it as a membrane peptide",
}


# ---------------------------------------------------------------- inputs

@dataclass
class Sources:
    """Everything the verdicts are built from. Loaders fill it; tests build it by hand."""
    curation: dict[str, dict] = field(default_factory=dict)  # curation.json "entries"
    mpstruc: dict[str, str] = field(default_factory=dict)  # PDB -> plain-English listing
    mpstruc_entries: set[str] = field(default_factory=set)  # parents and daughters (not related)
    pdbtm: set[str] = field(default_factory=set)
    memprotmd: set[str] = field(default_factory=set)
    opm: dict[str, str] = field(default_factory=dict)  # membranedb.opm_spanning labels
    ours: set[str] = field(default_factory=set)  # in candidates.json or curation.json
    summaries: dict[str, dict | None] = field(default_factory=dict)  # PDB -> {"title", "released"}
    scores: dict[str, dict | None] = field(default_factory=dict)  # PDB -> {"p", "why", "facts"}
    scorer: Callable[[list[str]], dict[str, dict | None]] | None = None

    def score(self, ids) -> None:
        """Model scores for the IDs that have none yet (None = the PDB does not know it)."""
        todo = sorted({p for p in ids if p not in self.scores})
        if todo and self.scorer:
            self.scores.update(self.scorer(todo))


def _read_json(path: Path, default=None):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def _write_json(path: Path, data, compact: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, separators=(",", ":")) if compact else json.dumps(data, indent=1), encoding="utf-8")
    os.replace(tmp, path)


def snapshot_from_export(path: Path) -> dict:
    """A minimal vercel_snapshot.json built straight from a Vercel export CSV (fallback only)."""
    day = (re.search(r"(\d{4}-\d{2}-\d{2})", path.name) or [None, dt.date.today().isoformat()])[1]
    entries = {}
    with path.open(newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            pdb = (r.get("Pdb Code") or "").strip().upper()
            if not pdb:
                continue
            bib = (r.get("Bibliography") or "") + " " + (r.get("PubMed Link") or "")
            pm = re.search(r"PMID:\s*(\d+)|term=(\d+)", bib)
            entries[pdb] = {
                "status": STATUS_FROM_EXPORT.get((r.get("Status") or "").strip().lower(), "pending"),
                "family": (r.get("Family") or "").strip() or None, "name": (r.get("Name") or "").strip() or None,
                "group": (pm.group(1) or pm.group(2)) if pm else None, "released": None, "tm": None,
                "memo": (r.get("Memo") or "").strip() or None, "seenAt": day}
    return {"source": f"export:{path.name}", "readAt": f"{day}T00:00:00Z", "entries": entries}


def load_snapshot(path: Path = SNAPSHOT) -> dict:
    snap = _read_json(path)
    if snap and snap.get("entries") is not None:
        return snap
    exports = sorted(EXPORTS.glob("*.csv"))
    if not exports:
        raise SystemExit(f"No Vercel snapshot: {path} is missing and {EXPORTS} has no export")
    print(f"{path.name} not found; reading {exports[-1].name} directly")
    return snapshot_from_export(exports[-1])


def mpstruc_listing(db) -> tuple[dict[str, str], set[str]]:
    """({PDB: 'Listed in the official mpstruc (...)'}, PDB codes of its entries)."""
    out: dict[str, str] = {}
    main: set[str] = set()
    for sg in db.subgroups():
        for e in sg.parents:
            for x in [e, *e.daughters]:
                code = x.pdbCode.upper()
                main.add(code)
                out[code] = f"Listed in the official mpstruc ({sg.name})"
                for r in x.related:
                    out.setdefault(r.pdb.upper(), f"Listed in the official mpstruc as related to {code} ({sg.name})")
    return out, main


def annotated_cached(start: str, end: str, kind: str, max_age_s: float = 7 * 86400) -> set[str]:
    """Entries released in [start, end] that the PDB tags with PDBTM or MemProtMD (cached a week)."""
    from .membrane import annotated
    cache = _read_json(ANNOT_CACHE, {})
    key = f"{kind}:{start}:{end}"
    hit = cache.get(key)
    if hit and time.time() - hit["at"] < max_age_s:
        return set(hit["ids"])
    ids = annotated(start, end, [kind])
    cache[key] = {"at": time.time(), "ids": sorted(ids)}
    _write_json(ANNOT_CACHE, cache)
    return ids


def model_scorer(db, cache_path: Path = SCORE_CACHE, chunk: int = 1000) -> Callable:
    """A scorer doing what membrane.score_structures does, with one Context for all
    batches, the resolved-segment facts kept for the reasons, and every batch saved
    to cache_path (keyed to the model file) so a long run resumes where it stopped."""
    from . import membrane

    state: dict = {}

    def score(ids: list[str]) -> dict[str, dict | None]:
        model = membrane.load_model()
        if not model:
            return {}
        key = hashlib.sha1(membrane.MODEL_PATH.read_bytes()).hexdigest()[:12]
        cache = _read_json(cache_path, {})
        if cache.get("model") != key:
            cache = {"model": key, "scores": {}}
        todo = [p for p in ids if p not in cache["scores"]]
        for i in range(0, len(todo), chunk):
            part = todo[i:i + chunk]
            if "ctx" not in state:
                state["ctx"] = membrane.Context(db, model.get("fam_stats"))
            ext, fam, det, cov, seg = membrane._fetch(part)
            for pdb in part:
                f = membrane.raw_facts(pdb, ext.get(pdb), fam.get(pdb), det.get(pdb), cov.get(pdb), seg)
                if not f:
                    cache["scores"][pdb] = None
                    continue
                x, why = membrane.features(f, state["ctx"], model["prior"])
                p = membrane.predict(x, model)
                cache["scores"][pdb] = {
                    "p": round(p, 3), "why": membrane.explain(x, why, f, p, model),
                    "facts": {k: f.get(k) for k in ("has_tm_protein", "has_barrel", "tm_in_construct", "tm_modelled",
                                                    "barrel_in_construct", "barrel_modelled", "tm_unmapped")}}
            _write_json(cache_path, cache)
            print(f"  scored {min(i + chunk, len(todo))}/{len(todo)}")
        return {p: cache["scores"].get(p) for p in ids}

    return score


def load_sources(db=None) -> Sources:
    from . import membranedb
    from .model import load
    from .report import source_xml
    db = db or load(source_xml())
    cur = (_read_json(SITE_DATA / "curation.json", {}) or {}).get("entries") or {}
    cands = (_read_json(SITE_DATA / "candidates.json", {}) or {}).get("candidates") or []
    listing, main = mpstruc_listing(db)
    src = Sources(curation={k.upper(): v for k, v in cur.items()}, mpstruc=listing, mpstruc_entries=main,
                  opm=membranedb.opm_spanning(membranedb.opm_detail()), ours={c["pdb"].upper() for c in cands} | {k.upper() for k in cur},
                  scorer=model_scorer(db))
    return src


def add_summaries(src: Sources, ids) -> None:
    from . import rcsb
    todo = sorted({p for p in ids if p not in src.summaries})
    if not todo:
        return
    for p, s in rcsb.entry_summaries(todo).items():
        src.summaries[p] = None if not s else {
            "title": (s.get("struct") or {}).get("title"),
            "released": ((s.get("rcsb_accession_info") or {}).get("initial_release_date") or "")[:10] or None}


# ---------------------------------------------------------------- verdicts

def segment_fact(facts: dict | None) -> str | None:
    """The resolved-segment facts in plain English."""
    if not facts:
        return None
    n, k = facts.get("tm_in_construct"), facts.get("tm_modelled")
    if facts.get("barrel_modelled"):
        return "the membrane-spanning beta-barrel is resolved"
    if k:
        return f"{k} membrane-crossing segment{'' if k == 1 else 's'} resolved" + (f" (of {n} in the construct)" if n and n > k else "")
    if n:
        return f"none of the {n} membrane-crossing segment{'' if n == 1 else 's'} in the construct is resolved"
    if facts.get("has_tm_protein") or facts.get("has_barrel"):
        return "the construct leaves out the membrane-crossing part of its proteins"
    if facts.get("has_tm_protein") is False and not facts.get("has_barrel"):
        return "UniProt marks no membrane-crossing segment in any of its proteins"
    return None


def model_text(s: dict) -> str:
    pct = round(s["p"] * 100)
    fact = segment_fact(s.get("facts"))
    why = [w for w in s.get("why") or [] if not (fact and "membrane-crossing segment" in w and "resolved" in w)]
    parts = ([fact] if fact else []) + why[:2 if fact else 3]
    return f"Our model: {pct}%" + (" - " + "; ".join(parts) if parts else "")


def membrane_part(facts: dict | None) -> bool | None:
    """Whether a membrane-crossing segment or barrel is resolved in the structure
    (None when the facts cannot tell: not scored, coverage unknown, or proteins
    UniProt calls transmembrane without mapping their segments)."""
    if not facts or facts.get("tm_modelled") is None:
        return None
    if facts.get("tm_modelled") or facts.get("barrel_modelled"):
        return True
    if facts.get("tm_unmapped"):
        return None
    return False


def mpstruc_doubt(pdb: str, src: Sources) -> str | None:
    """Why the official mpstruc listing of pdb is not trusted, or None.

    'no_membrane_protein': no protein in the structure at this code is a membrane
    protein at all (UniProt marks no membrane-crossing segment or barrel): a soluble
    part or partner of a membrane complex (F1, V1), a monotopic protein, or a typo in
    mpstruc's code (6IJ6 for 6I6J). 'fragment': its proteins are membrane proteins,
    but the structure leaves out every membrane-crossing part."""
    if pdb not in src.mpstruc:
        return None
    f = (src.scores.get(pdb) or {}).get("facts")
    if membrane_part(f) is not False:
        return None
    return "no_membrane_protein" if f.get("has_tm_protein") is False and not f.get("has_barrel") else "fragment"


MPSTRUC_DOUBT_TEXT = {
    "no_membrane_protein": "but no protein in the structure at this PDB code is a membrane protein (a soluble part "
                           "or partner, or a possible mpstruc ID typo), so mpstruc gives no vote",
    "fragment": "but mpstruc lists a fragment without the membrane part (no membrane-crossing segment or barrel "
                "is resolved), so mpstruc gives no vote",
}


def judge(pdb: str, src: Sources) -> tuple[str, list[dict], float, str]:
    """(verdict, evidence strongest first, confidence used for sorting, deciding source or 'none')."""
    ev: list[dict] = []
    verdict, conf, by = None, 0.0, "none"
    c = src.curation.get(pdb) or {}
    if c.get("status") == "approved":
        ev.append({"source": "curators", "text": "Our curators approved it for mpstruc"})
        verdict, conf, by = MEMBRANE, 3.0, "curators"
    elif c.get("status") == "rejected":
        reason = c.get("reason") or c.get("reasonCode") or "no reason given"
        if c.get("reasonCode") == "not_membrane":
            ev.append({"source": "curators", "text": f"Our curators rejected it: {reason}"})
            verdict, conf, by = NOT_MEMBRANE, 3.0, "curators"
        else:
            ev.append({"source": "curators", "text": f"Our curators set it aside: {reason} "
                                                     "(this says nothing about the membrane)"})
    s = src.scores.get(pdb)
    facts = (s or {}).get("facts")
    if pdb in src.mpstruc:
        doubt = mpstruc_doubt(pdb, src)
        ev.append({"source": "mpstruc", "text": src.mpstruc[pdb] + (f", {MPSTRUC_DOUBT_TEXT[doubt]}" if doubt else "")})
        if not verdict and not doubt:
            verdict, conf, by = MEMBRANE, 2.5, "mpstruc"
    dbs, mpm_only_weak = 0, False
    if pdb in src.pdbtm:
        ev.append({"source": "pdbtm", "text": "PDBTM lists it as transmembrane"})
        dbs += 1
    opm_tm = src.opm.get(pdb) == "transmembrane"
    if pdb in src.memprotmd:
        # MemProtMD simulates any structure in a membrane; on its own it is not enough when
        # nothing membrane-crossing is resolved and our model finds no sign of a membrane protein.
        weak = (not dbs and not opm_tm and s is not None and s["p"] < P_NOT and not membrane_part(facts))
        if weak:
            mpm_only_weak = True
            ev.append({"source": "memprotmd", "text": "MemProtMD lists it as a membrane protein, but nothing "
                                                      "membrane-crossing is resolved and our model finds no sign of one, "
                                                      "so this alone is not enough"})
        else:
            ev.append({"source": "memprotmd", "text": "MemProtMD lists it as a membrane protein"})
            dbs += 1
    if src.opm.get(pdb):
        ev.append({"source": "opm", "text": OPM_TEXT.get(src.opm[pdb], f"OPM lists it ({src.opm[pdb]})")})
        dbs += opm_tm
    if dbs and not verdict:
        verdict, conf, by = MEMBRANE, 2.0 + 0.1 * dbs, "databases"
    if s:
        ev.append({"source": "model", "text": model_text(s)})
        if not verdict:
            p = s["p"]
            verdict, conf = (MEMBRANE, p) if p >= P_MEMBRANE else (NOT_MEMBRANE, 1 - p) if p < P_NOT else (UNCERTAIN, 0.0)
            if verdict == NOT_MEMBRANE and mpm_only_weak:
                verdict, conf = UNCERTAIN, 0.0
            by = "model" if verdict != UNCERTAIN else "none"
    return verdict or UNCERTAIN, ev, conf, by


def decided_without_model(pdb: str, src: Sources) -> bool:
    """True when the verdict does not depend on our model's score or the resolved-segment
    facts that come with it (the official mpstruc listing and MemProtMD are checked against
    those facts, so they need a score)."""
    c = src.curation.get(pdb) or {}
    return (c.get("status") == "approved" or (c.get("status") == "rejected" and c.get("reasonCode") == "not_membrane")
            or pdb in src.pdbtm or src.opm.get(pdb) == "transmembrane")


def row(pdb: str, src: Sources, snap: dict) -> tuple[dict, float, str]:
    """(ROW, sorting confidence, deciding source)."""
    verdict, ev, conf, by = judge(pdb, src)
    v = snap.get(pdb) or {}
    s = src.summaries.get(pdb) or {}
    sc = src.scores.get(pdb)
    cur = src.curation.get(pdb) or {}
    return {"pdb": pdb, "title": s.get("title") or v.get("name") or cur.get("title"),
            "released": s.get("released") or v.get("released") or cur.get("released"),
            "verdict": verdict, "evidence": ev, "ourScore": sc["p"] if sc else None,
            "vercelStatus": v.get("status") if pdb in snap else None,
            "vercelFamily": v.get("family") if pdb in snap else None,
            "inMpstruc": pdb in src.mpstruc, "inOurList": pdb in src.ours}, conf, by


def _sorted(rows: list[tuple]) -> list[dict]:
    """Most confident first: deciding source, then our model's score (highest first for
    membrane rows, lowest first for the others), then PDB ID."""
    def key(t):
        r, conf = t[0], t[1]
        sc = r["ourScore"]
        tie = -(sc or 0) if r["verdict"] == MEMBRANE else (sc if sc is not None else 1)
        return -conf, tie, r["pdb"]
    return [t[0] for t in sorted(rows, key=key)]


# ---------------------------------------------------------------- comparisons

def release_wednesday(day: dt.date) -> dt.date:
    """The latest PDB release day (Wednesday) on or before day."""
    return day - dt.timedelta(days=(day.weekday() - 2) % 7)


def week_window(week: str) -> tuple[str, str]:
    end = dt.date.fromisoformat(week)
    return (end - dt.timedelta(days=6)).isoformat(), week


def compare_week(week: str, ids: list[str], src: Sources, snap: dict) -> dict:
    """One release week: ids are the PDB entries released in it."""
    start, end = week_window(week)
    ids = sorted({p.upper() for p in ids})
    need = [p for p in ids if not decided_without_model(p, src)]
    src.score(need)
    listed = set(ids) & snap.keys() | {p for p, v in snap.items() if v.get("released") and start <= v["released"] <= end}
    missed, notmem, agree, verdicts = [], [], 0, {MEMBRANE: 0, NOT_MEMBRANE: 0, UNCERTAIN: 0}
    for p in ids:
        r = row(p, src, snap)
        verdicts[r[0]["verdict"]] += 1
        if p in snap:
            if r[0]["verdict"] == MEMBRANE:
                agree += 1
            elif r[0]["verdict"] == NOT_MEMBRANE:
                notmem.append(r)
        elif r[0]["verdict"] == MEMBRANE:
            missed.append(r)
    return {"week": week, "from": start, "released": len(ids), "verdicts": verdicts, "vercelListed": len(listed),
            "missed": _sorted(missed), "notMembrane": _sorted(notmem), "agree": agree}


def compare_whole(snap: dict, src: Sources, span_ids: list[str], recent_ids: list[str]) -> dict:
    """Vercel's whole snapshot. span_ids: PDB entries released within Vercel's span;
    recent_ids: those of the last RECENT_WEEKS weeks of it (our model is used there)."""
    vids = sorted(snap)
    src.score([p for p in vids if not decided_without_model(p, src)])
    notmem, by_source, verdicts = [], {}, {MEMBRANE: 0, NOT_MEMBRANE: 0, UNCERTAIN: 0}
    said = {}
    for p in vids:
        r, c, by = row(p, src, snap)
        verdicts[r["verdict"]] += 1
        by_source[by] = by_source.get(by, 0) + 1
        if r["verdict"] == NOT_MEMBRANE:
            notmem.append((r, c))
            said[r["vercelStatus"]] = said.get(r["vercelStatus"], 0) + 1
    span = set(p.upper() for p in span_ids)
    strong = {p for p in span if p in src.mpstruc or p in src.pdbtm or p in src.memprotmd
              or src.opm.get(p) == "transmembrane"}
    recent = {p.upper() for p in recent_ids} - snap.keys() - strong
    absent = sorted(set(src.mpstruc) - snap.keys())
    # Strong-source candidates the mpstruc and MemProtMD checks need a score for, and the
    # mpstruc codes absent from Vercel (checked for typos and fragments).
    src.score([p for p in recent | (strong - snap.keys()) | set(absent) if not decided_without_model(p, src)])
    pool = (strong - snap.keys()) | {p for p in recent if (src.scores.get(p) or {}).get("p", 0) >= P_MEMBRANE}
    missed = [r for r in (row(p, src, snap) for p in sorted(pool)) if r[0]["verdict"] == MEMBRANE]
    missed_by: dict[str, int] = {}  # the source that decided each missed row (databases split by the first one)
    for r, _, by in missed:
        k = next((e["source"] for e in r["evidence"] if e["source"] in ("pdbtm", "memprotmd", "opm")), by) \
            if by == "databases" else by
        missed_by[k] = missed_by.get(k, 0) + 1
    already = sorted(p for p, v in snap.items() if v.get("status") == "already_in" and p not in src.mpstruc)
    # A code mpstruc has that Vercel lacks, and a code Vercel marks 'Already in' that mpstruc
    # lacks, with the same characters in two places swapped (7JOI / 7OJI), are most likely
    # one structure with a typo on one side.
    # (Only 'Already in' codes: pairing with any Vercel code finds chance look-alikes.)
    doubt = {p: mpstruc_doubt(p, src) for p in absent}
    pairs = typo_pairs(absent, already)
    paired_m, paired_v = {a for a, _ in pairs}, {b for _, b in pairs}
    doubtful = sorted(p for p in absent if doubt[p] and p not in paired_m)
    already = [p for p in already if p not in paired_v]
    absent = [p for p in absent if p not in paired_m and not doubt[p]]
    absent_main = [p for p in absent if p in src.mpstruc_entries]
    summary = {
        "snapshotEntries": len(vids), "verdicts": verdicts, "decidedBy": by_source,
        "scored": sum(1 for p in vids if src.scores.get(p)),
        "unknownToPdb": sum(1 for p in vids if p in src.scores and src.scores[p] is None
                            and not decided_without_model(p, src)),
        "notMembrane": len(notmem), "notMembraneByVercelStatus": said,
        "notMembraneVercelSaidNo": said.get("no", 0),
        "missed": len(missed), "missedStrongSources": sum(1 for r, *_ in missed if r["pdb"] in strong),
        "missedModelRecent": sum(1 for r, *_ in missed if r["pdb"] not in strong),
        "missedBySource": missed_by,
        "spanReleased": len(span), "recentScored": len(recent),
        "alreadyInNotInMpstruc": {"count": len(already), "pdbs": already},
        "mpstrucNotInVercel": {"count": len(absent), "entries": len(absent_main),
                               "relatedOnly": len(absent) - len(absent_main),
                               "inVercelSpan": sum(1 for p in absent if p in span), "pdbs": absent},
        "possibleMpstrucTypos": {"count": len(pairs),
                                 "pairs": [{"mpstruc": a, "vercel": b, "vercelStatus": snap[b].get("status"),
                                            "mpstrucDoubt": doubt.get(a)} for a, b in pairs]},
        "mpstrucWithoutMembranePart": {
            "count": len(doubtful), "pdbs": doubtful,
            "noMembraneProtein": [p for p in doubtful if doubt[p] == "no_membrane_protein"],
            "fragments": [p for p in doubtful if doubt[p] == "fragment"]},
    }
    return {"summary": summary, "missed": _sorted(missed), "notMembrane": _sorted(notmem)}


def swapped(a: str, b: str) -> bool:
    """True when b is a with two of its characters swapped (6IJ6 / 6I6J, 6SMX / 6XMS)."""
    if len(a) != len(b) or a == b:
        return False
    diff = [i for i in range(len(a)) if a[i] != b[i]]
    return len(diff) == 2 and a[diff[0]] == b[diff[1]] and a[diff[1]] == b[diff[0]]


def typo_pairs(mpstruc_only: list[str], vercel_only: list[str]) -> list[tuple[str, str]]:
    """(mpstruc code, Vercel code) pairs that differ by two swapped characters."""
    by_key: dict[str, list[str]] = {}
    for b in vercel_only:
        by_key.setdefault("".join(sorted(b)), []).append(b)
    return sorted((a, b) for a in mpstruc_only for b in by_key.get("".join(sorted(a)), []) if swapped(a, b))


def snapshot_counts(snap_doc: dict) -> dict:
    e = snap_doc.get("entries") or {}
    by_status, by_family = {}, {}
    for v in e.values():
        by_status[v.get("status") or "pending"] = by_status.get(v.get("status") or "pending", 0) + 1
        fam = v.get("family") or "(blank)"
        by_family[fam] = by_family.get(fam, 0) + 1
    rel = sorted(v["released"] for v in e.values() if v.get("released"))
    out = {"source": snap_doc.get("source"), "readAt": snap_doc.get("readAt"), "entries": len(e),
           "byStatus": by_status, "byFamily": by_family}
    if rel:
        out["releasedSpan"] = [rel[0], rel[-1]]
    return out


def merge(old: dict | None, snapshot: dict, week: dict | None = None, whole: dict | None = None) -> dict:
    old = old or {}
    weeks = [w for w in old.get("weeks") or [] if not week or w.get("week") != week["week"]]
    if week:
        weeks.append(week)
    weeks.sort(key=lambda w: w["week"], reverse=True)
    out = {"generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "snapshot": snapshot,
           "weeks": weeks}
    if whole or old.get("wholeDatabase"):
        out["wholeDatabase"] = whole or old["wholeDatabase"]
    return out


# ---------------------------------------------------------------- CLI

def span_of(snap: dict, src: Sources) -> tuple[str, str]:
    add_summaries(src, list(snap))
    days = sorted(d for d in ((src.summaries.get(p) or {}).get("released") or (snap[p].get("released"))
                              for p in snap) if d)
    if not days:
        raise SystemExit("No release dates for Vercel's snapshot")
    return days[0], days[-1]


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--week", help="release Wednesday YYYY-MM-DD (default: the latest one)")
    ap.add_argument("--whole", action="store_true", help="compare Vercel's whole snapshot")
    ap.add_argument("--snapshot", type=Path, default=SNAPSHOT)
    ap.add_argument("--out", type=Path, default=OUTPUT)
    args = ap.parse_args(argv)
    from .automation import released_between

    snap_doc = load_snapshot(args.snapshot)
    snap = {k.upper(): v for k, v in (snap_doc.get("entries") or {}).items()}
    src = load_sources()
    week = whole = None
    if args.week or not args.whole:
        wk = args.week or release_wednesday(dt.datetime.now(dt.timezone.utc).date()).isoformat()
        if dt.date.fromisoformat(wk).weekday() != 2:
            raise SystemExit(f"{wk} is not a Wednesday (the PDB releases on Wednesdays)")
        start, end = week_window(wk)
        ids = released_between(start, end)
        src.pdbtm |= annotated_cached(start, end, "PDBTM")
        src.memprotmd |= annotated_cached(start, end, "MemProtMD")
        add_summaries(src, ids)
        week = compare_week(wk, ids, src, snap)
        print(f"week {wk}: {week['released']} released, Vercel lists {week['vercelListed']}, "
              f"{len(week['missed'])} missed, {len(week['notMembrane'])} not membrane, {week['agree']} agree")
    if args.whole:
        lo, hi = span_of(snap, src)
        print(f"Vercel's snapshot: {len(snap)} entries released {lo} .. {hi}")
        span_ids = released_between(lo, hi)
        recent_lo = max(lo, (dt.date.fromisoformat(hi) - dt.timedelta(weeks=RECENT_WEEKS)).isoformat())
        recent_ids = released_between(recent_lo, hi)
        src.pdbtm |= annotated_cached("1970-01-01", hi, "PDBTM")
        src.memprotmd |= annotated_cached("1970-01-01", hi, "MemProtMD")
        whole = compare_whole(snap, src, span_ids, recent_ids)
        rows = whole["missed"] + whole["notMembrane"]
        add_summaries(src, [r["pdb"] for r in rows])
        for r in rows:  # titles and dates for the rows only
            s = src.summaries.get(r["pdb"]) or {}
            r["title"] = s.get("title") or r["title"]
            r["released"] = s.get("released") or r["released"]
        whole["summary"]["span"] = [lo, hi]
        whole["summary"]["recentSpan"] = [recent_lo, hi]
        print(json.dumps({k: v for k, v in whole["summary"].items()
                          if k not in ("alreadyInNotInMpstruc", "mpstrucNotInVercel")}))
        print(f"'Already in' at Vercel but not in mpstruc: {whole['summary']['alreadyInNotInMpstruc']['count']}; "
              f"mpstruc codes absent from Vercel: {whole['summary']['mpstrucNotInVercel']['count']}")
    out = merge(_read_json(args.out), snapshot_counts(snap_doc), week, whole)
    _write_json(args.out, out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
