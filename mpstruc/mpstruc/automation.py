"""Find newly released membrane-protein structures and propose mpstruc placements.

    python -m mpstruc.automation run [--since YYYY-MM-DD] [--until YYYY-MM-DD]
    python -m mpstruc.automation backtest

`run` writes ../data/candidates.json for the site's Automation tab.
`backtest` replays past periods against mpstruc's own decisions, fits the
likelihood model and cut-off on the first, reports on the second, and
writes model.json, reports/backtest.{md,json} and ../data/backtest.json.

Every decision carries its evidence and a confidence taken from the
backtest. The rules come from reports/placement_analysis.md:
  related   same primary citation as an existing entry, same protein   (93%)
  daughter  same protein(s) as an earlier entry, different paper       (94%)
  partial   shares some proteins only (complexes, fusions)             (67%: review)
  parent    protein not yet in mpstruc                                 (90%)
RCSB's copy of the mpstruc classification is never used as evidence.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from . import rcsb
from .model import Database, load
from .placement import MEMBRANE_SOURCES, entity_summary
from .report import ROOT, source_xml
from .text import to_text

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
SITE_DATA = ROOT.parent / "data"
# Fusion partners and crystallisation chaperones that make a UniProt set look
# different although the membrane protein is the same.
FUSION_UNIPROT = {
    "P00720",  # T4 lysozyme
    "P0ABE7",  # BRIL (apocytochrome b562RIL)
    "P0AEX9",  # maltose-binding protein
    "P42212",  # GFP
    "P0AA25",  # thioredoxin
    "Q9ZEU2",  # PGS (Pyrococcus abyssi glycogen synthase)
    "P61626",  # human lysozyme
    "P63165",  # SUMO1
    "P0ACJ8",  # rubredoxin
    "P02185",  # myoglobin
}


# ---------------------------------------------------------------- PDB search

def released_between(start: str, end: str) -> list[str]:
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_accession_info.initial_release_date", "operator": "range",
            "value": {"from": start, "to": end, "include_lower": True, "include_upper": True}}},
         "return_type": "entry", "request_options": {"return_all_hits": True}}
    req = urllib.request.Request(SEARCH, data=json.dumps(q).encode(),
                                 headers={"Content-Type": "application/json", "User-Agent": "mpstruc-curation-tool"})
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            break
        except (urllib.error.URLError, TimeoutError, ConnectionError):
            if attempt == 5:
                raise
            import time
            time.sleep(2 ** (attempt + 1))
    if not body:  # 204: no hits
        return []
    return sorted(x["identifier"] for x in json.loads(body).get("result_set", []))


# ---------------------------------------------------------------- knowledge base

@dataclass
class Known:
    pdb: str
    key: str
    kind: str  # parent | daughter | related
    sg: str  # subgroup key
    parent_key: str  # owning parent record key
    owner_key: str  # entry the related ID is listed under (== key for parents/daughters)
    released: str
    s: dict | None


@dataclass
class KB:
    """What mpstruc knew at a given date."""
    entries: list[Known]
    by_pdb: dict[str, list[Known]] = field(default_factory=dict)
    by_key: dict[str, Known] = field(default_factory=dict)
    by_uniprot: dict[str, list[Known]] = field(default_factory=dict)
    by_citation: dict[str, list[Known]] = field(default_factory=dict)
    by_authors: dict[str, list[Known]] = field(default_factory=dict)
    feat_sg: dict[str, Counter] = field(default_factory=dict)  # feature -> subgroup counts
    sg_count: Counter = field(default_factory=Counter)
    sg_names: dict[str, str] = field(default_factory=dict)
    sg_group: dict[str, str] = field(default_factory=dict)


def _features(s: dict | None) -> set[str]:
    if not s:
        return set()
    return {f"pfam:{x}" for x in s["pfam"]} | {f"ipr:{x}" for x in s.get("interpro", ())} | \
           {f"opm:{x}" for x in s.get("opm", ())} | {f"up:{x}" for x in s["core"]} | \
           {f"uipr:{x}" for x in s.get("uipr", ())}


TM_KEYWORDS = {"KW-0812", "KW-1133", "KW-1134"}  # Transmembrane, TM helix, TM beta strand


def _authors_key(authors: list[str] | None) -> str | None:
    """First and last depositing author, e.g. 'auth:ahmed|lee': a stand-in for
    'same paper' when the PDB has not linked the paper yet."""
    names = [a.split(",")[0].strip().lower() for a in authors or [] if a.strip()]
    return f"auth:{names[0]}|{names[-1]}" if names else None


def summarise(entry: dict | None, extra: dict | None = None, fam: dict | None = None) -> dict | None:
    """What the automation knows about a structure.

    `core` is the structure's membrane protein(s): the UniProt accessions
    whose UniProt entry carries a transmembrane keyword, or, if none does,
    all of them; fusion partners are dropped. This uses UniProt rather than
    the PDB's OPM/PDBTM/MemProtMD tags, which new releases no longer get, so
    old and new structures are described the same way.
    """
    s = entity_summary(entry)
    if s is None:
        return None
    anns = [a for e in (entry.get("polymer_entities") or [])
            for a in e.get("rcsb_polymer_entity_annotation") or []]
    s["interpro"] = {a["annotation_id"] for a in anns if a.get("type") == "InterPro"}
    s["opm"] = {".".join(a["annotation_id"].split(".")[:3]) for a in anns if a.get("type") == "OPM"}
    kw = (extra or {}).get("kw") or {}
    tm = {u for u in s["all_uniprot"] if TM_KEYWORDS & set(kw.get(u, ()))}
    s["core"] = (tm or s["all_uniprot"]) - FUSION_UNIPROT
    s["uniprot"] = s["core"] | (s["all_uniprot"] & FUSION_UNIPROT)
    s["keywords"] = {k for ks in kw.values() for k in ks}
    ipr = (fam or {}).get("ipr") or {}
    s["uipr"] = {x for u in s["core"] for x in ipr.get(u, ())}
    s["authors_key"] = _authors_key((fam or {}).get("authors"))
    rel = s.get("released") or ""
    late = [d for d in (extra or {}).get("db_ref_updates", []) if rel and d <= _plus_days(rel, 365)]
    s["citation_late"] = bool(late)
    s["title"] = None
    return s


def _plus_days(day: str, n: int) -> str:
    return (dt.date.fromisoformat(day) + dt.timedelta(days=n)).isoformat()


def as_at_release(s: dict) -> dict:
    """The structure as the weekly job would have seen it in its release week:
    no paper if the PDB linked it later, and none of the annotations the PDB
    adds afterwards (OPM/PDBTM/MemProtMD, Pfam, InterPro)."""
    s = dict(s)
    if s.get("citation_late"):
        s["citation"], s["doi"], s["pmid"] = set(), None, None
    s["sources"], s["opm"], s["pfam"], s["interpro"] = set(), set(), set(), set()
    return s


def build_kb(db: Database, S: dict[str, dict | None], cutoff: str | None = None) -> KB:
    kb = KB(entries=[])
    for g in db.groups:
        for sg in g.subgroups:
            kb.sg_names[sg.key] = to_text(sg.name) or ""
            kb.sg_group[sg.key] = g.topology
            for p in sg.parents:
                for e in (p, *p.daughters):
                    items = [(e.pdbCode, e.key, e.kind, e.key)] + [(r.pdb, r.key, "related", e.key) for r in e.related]
                    for pdb, key, kind, owner in items:
                        s = S.get(pdb.upper())
                        rel = (s or {}).get("released") or ""
                        if cutoff and (not rel or rel >= cutoff):
                            continue
                        kb.entries.append(Known(pdb.upper(), key, kind, sg.key, p.key, owner, rel, s))
    for k in kb.entries:
        kb.by_pdb.setdefault(k.pdb, []).append(k)
        kb.by_key[k.key] = k
        if not k.s:
            continue
        for u in k.s["core"]:
            kb.by_uniprot.setdefault(u, []).append(k)
        for c in k.s["citation"]:
            kb.by_citation.setdefault(c, []).append(k)
        if k.s.get("authors_key"):
            kb.by_authors.setdefault(k.s["authors_key"], []).append(k)
        for f in _features(k.s):
            kb.feat_sg.setdefault(f, Counter())[k.sg] += 1
        kb.sg_count[k.sg] += 1
    return kb


# ---------------------------------------------------------------- classification

def subgroup_scores(kb: KB, s: dict) -> list[tuple[str, float, list[str]]]:
    """Rank subgroups for a structure by its UniProt, Pfam, InterPro and OPM family.

    Each feature votes for the subgroups where it occurs, weighted by how
    specific it is (a feature seen in one subgroup counts more than one
    seen in twenty). Returns [(subgroup, share of the vote, evidence)].
    """
    weights = {"up": 4.0, "opm": 2.0, "pfam": 1.5, "uipr": 1.5, "ipr": 1.0}
    votes: Counter = Counter()
    why: dict[str, list[str]] = defaultdict(list)
    for f in _features(s):
        dist = kb.feat_sg.get(f)
        if not dist:
            continue
        total = sum(dist.values())
        spec = 1.0 / math.log2(1 + len(dist))
        for sg, n in dist.items():
            votes[sg] += weights[f.split(":")[0]] * spec * n / total
            if len(why[sg]) < 4:
                why[sg].append(f)
    tot = sum(votes.values())
    return [(sg, v / tot, why[sg]) for sg, v in votes.most_common(3)] if tot else []


def _owner(kb: KB, k: Known) -> Known:
    """The parent/daughter record a known entry hangs under (itself unless related)."""
    return kb.by_key.get(k.owner_key, k)


@dataclass
class Proposal:
    pdb: str
    placement: str  # parent | daughter | related | review
    subgroup: str | None
    alternatives: list
    owner_pdb: str | None  # for daughter: the parent; for related: the entry it is listed under
    owner_key: str | None
    confidence: float
    rule: str
    evidence: list[str]
    batch_lead: str | None = None


def classify(kb: KB, pdb: str, s: dict, batch_leads: dict[str, str], conf: dict[str, float]) -> Proposal:
    ev: list[str] = []
    core = s["core"]
    sgs = subgroup_scores(kb, s)
    alts = [{"subgroup": sg, "name": kb.sg_names[sg], "score": round(v, 3), "evidence": w} for sg, v, w in sgs]

    # 1. Related: same paper as an existing entry, same protein. A paper's
    #    later structures join the entry its earlier ones are listed under.
    same_paper = [k for c in s["citation"] for k in kb.by_citation.get(c, [])]
    same_paper = [k for k in same_paper if k.s and (k.s["core"] & core)]
    if not s["citation"] and s.get("authors_key") and s.get("released"):
        # Paper not linked yet: same depositing authors and protein, released
        # within six months, is most likely the same study.
        soon = _plus_days(s["released"], -183)
        same_paper = [k for k in kb.by_authors.get(s["authors_key"], [])
                      if k.s and (k.s["core"] & core) and k.released >= soon]
    if same_paper:
        owners = Counter(_owner(kb, k).key for k in same_paper)
        o = kb.by_key[sorted(owners, key=lambda key: (-owners[key], key))[0]]
        ev.append(f"Same primary citation and protein as {', '.join(sorted({k.pdb for k in same_paper})[:5])}, "
                  f"listed under {o.pdb} ({o.kind}, {o.key})")
        return Proposal(pdb, "related", o.sg, alts, o.pdb, o.key, conf.get("related", 0.9), "same-paper", ev)

    # 1b. Same paper and a shared protein with other new structures this week:
    #     related to the lowest ID among them (the batch lead).
    for c in paper_keys(s):
        members = [p for p, cr in batch_leads.get(c, []) if cr & core]
        lead = min(members) if len(members) > 1 else None
        if lead and lead != pdb:
            ev.append(f"Same primary citation as {lead}, released in the same batch")
            return Proposal(pdb, "related", None, alts, lead, None, conf.get("related-batch", 0.85),
                            "same-paper-batch", ev, batch_lead=lead)

    # 2. Same protein(s) as an earlier entry: daughter of that entry's parent.
    prior = [k for u in core for k in kb.by_uniprot.get(u, [])]
    exact = [k for k in prior if k.s and k.s["core"] == core]
    if exact:
        # The parent that already holds the most structures of this protein (then the earliest).
        votes = Counter(k.parent_key for k in exact)
        parent_key = sorted(votes, key=lambda key: (-votes[key], key))[0]
        k = next(k for k in exact if k.parent_key == parent_key)
        ev.append(f"Same protein(s) {', '.join(sorted(core))} as {votes[parent_key]} structure(s) under "
                  f"{_parent_pdb(kb, parent_key)} ({parent_key})")
        return Proposal(pdb, "daughter", k.sg, alts, _parent_pdb(kb, parent_key), parent_key,
                        conf.get("daughter", 0.9), "same-protein", ev)

    if prior:
        k = Counter((x.parent_key, x.sg) for x in prior).most_common(1)[0][0]
        shared = sorted(core & {u for x in prior for u in x.s["core"]})
        ev.append(f"Shares {', '.join(shared)} with {len(prior)} known entries (complex or fusion differs)")
        return Proposal(pdb, "review", k[1], alts, _parent_pdb(kb, k[0]), k[0], conf.get("partial", 0.6),
                        "partial-protein", ev)

    # 3. New protein: new parent in the best-scoring subgroup.
    if sgs:
        ev.append(f"No known entry shares its protein; subfamily from {', '.join(sgs[0][2])}")
        return Proposal(pdb, "parent", sgs[0][0], alts, None, None,
                        conf.get("parent", 0.8) * min(1.0, sgs[0][1] / 0.5), "new-protein", ev)
    ev.append("No UniProt, Pfam, InterPro or OPM family links it to any subfamily")
    return Proposal(pdb, "review", None, alts, None, None, 0.0, "no-evidence", ev)


def _parent_pdb(kb: KB, parent_key: str) -> str | None:
    k = kb.by_key.get(parent_key)
    return k.pdb if k else None


# OPM, PDBTM and MemProtMD tags are not used: the PDB stopped assigning them
# to new releases during 2025, so they never fire in weekly use.
SIGNALS = {
    "kw_tm": "UniProt: transmembrane protein",
    "kw_membrane": "UniProt: membrane protein",
    "same_protein": "Same protein as an mpstruc entry",
    "family": "Protein family seen in mpstruc",
    "same_paper": "Same paper as an mpstruc entry",
    "batch_paper": "Same paper as another candidate this week",
}


TITLE_STOP = {"of", "the", "in", "with", "and", "a", "an", "from", "to", "by", "for", "on", "at", "bound",
              "structure", "structures", "complex", "state", "form", "cryo-em", "cryoem", "crystal"}


def title_tokens(title: str | None) -> set[str]:
    import re
    words = re.findall(r"[a-z][a-z0-9-]{2,}", (title or "").lower())
    return {f"t:{w}" for w in words if w not in TITLE_STOP}


def signals(kb: KB, s: dict, batch_papers: set[str] = frozenset()) -> dict[str, float]:
    """Named membrane/mpstruc signals, plus method and title-word features for the model."""
    kw = s.get("keywords") or set()
    fam = [f for f in _features(s) if f.startswith(("pfam:", "opm:", "ipr:")) and f in kb.feat_sg]
    method = (s.get("method") or "").upper()
    extra = {
        "m:em": float("ELECTRON MICROSCOPY" in method),
        "m:xray": float("X-RAY" in method),
        "m:nmr": float("NMR" in method),
        "n_entities>3": float(s.get("n_protein_entities", 0) > 3),
        **{t: 1.0 for t in title_tokens(s.get("title"))},
    }
    return extra | {
        "kw_tm": float(bool(kw & TM_KEYWORDS)),
        "kw_membrane": float("KW-0472" in kw),
        "same_protein": float(any(u in kb.by_uniprot for u in s["core"])),
        "family": float(bool(fam)),
        "same_paper": float(any(c in kb.by_citation for c in s["citation"]) or
                            (not s["citation"] and bool(s.get("authors_key")) and
                             any(k.s and k.s["core"] & s["core"] for k in kb.by_authors.get(s["authors_key"], [])))),
        "batch_paper": float(bool(paper_keys(s) & batch_papers)),
    }


def in_pool(sig: dict[str, float]) -> bool:
    """Worth scoring at all: some membrane or mpstruc signal fired."""
    return any(sig.get(k) for k in ("kw_tm", "same_protein", "same_paper")) or \
        bool(sig.get("kw_membrane") and sig.get("family"))


def score(sig: dict[str, float], model: dict | None) -> float:
    """Probability that mpstruc lists the structure, from the logistic model."""
    if not model:
        return 0.5
    z = model["bias"] + sum(model["weights"].get(k, 0.0) * v for k, v in sig.items())
    return 1 / (1 + math.exp(-z))


def fit_model(X: list[dict[str, float]], y: list[int], l2: float = 6.0, steps: int = 4000,
              min_count: int = 15) -> dict:
    """Weighted logistic regression. Title words seen fewer than min_count times are ignored."""
    import numpy as np
    counts = Counter(k for x in X for k, v in x.items() if v)
    keys = sorted(k for k in counts if not k.startswith("t:") or counts[k] >= min_count)
    A = np.array([[x.get(k, 0.0) for k in keys] for x in X], dtype=float)
    t = np.array(y, dtype=float)
    w, b = np.zeros(len(keys)), 0.0
    pos = max(t.mean(), 1e-6)
    sw = np.where(t == 1, 0.5 / pos, 0.5 / (1 - pos))  # balance the rare positives
    for _ in range(steps):
        p = 1 / (1 + np.exp(-(A @ w + b)))
        g = sw * (p - t)
        w -= 0.5 * ((A.T @ g) / len(t) + l2 * w / len(t))
        b -= 0.5 * g.mean()
    return {"bias": round(float(b), 4),
            "weights": {k: round(float(v), 4) for k, v in zip(keys, w) if abs(v) >= 0.01}}


def paper_keys(s: dict) -> set[str]:
    """The paper a structure belongs to: its citation, or, until the PDB links
    one, its first and last depositing author."""
    return set(s["citation"]) or ({s["authors_key"]} if s.get("authors_key") else set())


def batch_leads(cands: dict[str, dict]) -> dict[str, list[tuple[str, set]]]:
    """New structures grouped by paper: {citation: [(pdb, proteins), ...]}.

    A paper's structures that share a protein are listed under the lowest
    PDB ID among them: in the snapshot the entry a paper's related structures
    hang under is the lowest ID 70% of the time, the best resolution 52%.
    """
    by_cit: dict[str, list[tuple[str, set]]] = defaultdict(list)
    for pdb, s in cands.items():
        for c in paper_keys(s):
            by_cit[c].append((pdb, s["core"]))
    return {c: v for c, v in by_cit.items() if len(v) > 1}


def _reasons(sig: dict[str, float], model: dict | None, n: int = 4) -> list[list]:
    """The features that moved the likelihood most, for display: [[label, weight], ...]."""
    if not model:
        return []
    contrib = [(k, model["weights"].get(k, 0.0)) for k, v in sig.items() if v and model["weights"].get(k)]
    contrib.sort(key=lambda kv: -abs(kv[1]))
    label = lambda k: SIGNALS.get(k) or {"m:em": "Cryo-EM", "m:xray": "X-ray", "m:nmr": "NMR",
                                         "n_entities>3": "More than 3 protein chains"}.get(k) or f"Title word '{k[2:]}'"
    return [[label(k), round(w, 2)] for k, w in contrib[:n]]


def propose(kb: KB, ids: list[str], conf: dict[str, float], model: dict | None = None,
            release_time: bool = False) -> tuple[list[dict], dict]:
    """Proposals for the given releases. release_time=True shows each structure
    as it looked in its release week (see as_at_release), for the backtest."""
    ents = rcsb.entry_entities(ids)
    summ = rcsb.entry_summaries(ids)
    extras = rcsb.entry_extras(ids)
    fams = rcsb.entry_families(ids)
    S = {}
    for pdb in ids:
        s = summarise(ents.get(pdb), extras.get(pdb), fams.get(pdb))
        if s is None:
            continue
        info = summ.get(pdb) or {}
        res = (info.get("rcsb_entry_info") or {}).get("resolution_combined") or []
        s["res"] = min(res) if res else None
        s["title"] = (info.get("struct") or {}).get("title")
        s["method"] = ", ".join(m.get("method", "") for m in info.get("exptl") or [])
        raw = ents.get(pdb) or {}
        cit = raw.get("rcsb_primary_citation") or {}
        s["doi"] = cit.get("pdbx_database_id_DOI")
        pm = cit.get("pdbx_database_id_PubMed")
        s["pmid"] = str(pm) if pm and int(pm) > 0 else None
        s["organisms"] = sorted({o.get("ncbi_scientific_name") for e in raw.get("polymer_entities") or []
                                 for o in e.get("rcsb_entity_source_organism") or [] if o.get("ncbi_scientific_name")})
        s["entities"] = [((e.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "")
                         for e in raw.get("polymer_entities") or []][:8]
        if release_time:
            s = as_at_release(s)
        S[pdb] = s
    stats = Counter()
    pool = {}
    for pdb, s in S.items():
        stats["released"] += 1
        if pdb in kb.by_pdb:
            stats["already listed"] += 1
            continue
        pool[pdb] = s
    # Papers with several new structures this batch (only among membrane-flagged ones).
    first = {pdb: signals(kb, s) for pdb, s in pool.items()}
    papers = Counter(c for pdb, s in pool.items() if in_pool(first[pdb]) for c in paper_keys(s))
    batch_papers = {c for c, n in papers.items() if n > 1}
    cands, sigs, scores = {}, {}, {}
    for pdb, s in pool.items():
        sig = signals(kb, s, batch_papers)
        if in_pool(sig):
            cands[pdb], sigs[pdb], scores[pdb] = s, sig, score(sig, model)
    leads = batch_leads(cands)
    out = []
    for pdb, s in sorted(cands.items()):
        p = classify(kb, pdb, s, leads, conf)
        if p.batch_lead and p.batch_lead in cands and p.batch_lead != pdb:
            lead = classify(kb, p.batch_lead, cands[p.batch_lead], {}, conf)
            p.subgroup = lead.subgroup
        out.append({
            "pdb": pdb, "title": s["title"], "method": s["method"], "res": s["res"], "released": s["released"],
            "uniprot": sorted(s["uniprot"]), "taxa": sorted(s["taxa"]), "organisms": s.get("organisms", []),
            "entities": s.get("entities", []), "doi": s.get("doi"), "pmid": s.get("pmid"),
            "signals": [SIGNALS[k] for k, v in sigs[pdb].items() if v and k in SIGNALS],
            "sig": sigs[pdb], "reasons": _reasons(sigs[pdb], model),
            "likelihood": round(scores[pdb], 3),
            "placement": p.placement, "subgroup": p.subgroup, "subgroupName": kb.sg_names.get(p.subgroup or ""),
            "group": kb.sg_group.get(p.subgroup or ""), "ownerPdb": p.owner_pdb, "ownerKey": p.owner_key,
            "confidence": round(p.confidence, 2), "rule": p.rule, "evidence": p.evidence,
            "alternatives": p.alternatives,
        })
    stats["candidates"] = len(out)
    return out, dict(stats)


# ---------------------------------------------------------------- backtest

def _truth(kb_full: KB, pdb: str) -> list[Known]:
    return kb_full.by_pdb.get(pdb, [])


def _full_summaries(db: Database) -> dict[str, dict | None]:
    ids = sorted({e.pdbCode.upper() for e in db.entries()} | {r.pdb.upper() for e in db.entries() for r in e.related})
    ents, extras, fams = rcsb.entry_entities(ids), rcsb.entry_extras(ids), rcsb.entry_families(ids)
    return {k: summarise(v, extras.get(k), fams.get(k)) for k, v in ents.items()}


def _weeks(start: str, end: str) -> list[tuple[str, str]]:
    d0, d1 = dt.date.fromisoformat(start), dt.date.fromisoformat(end)
    out = []
    while d0 <= d1:
        out.append((d0.isoformat(), min(d0 + dt.timedelta(days=6), d1).isoformat()))
        d0 += dt.timedelta(days=7)
    return out


def evaluate(db: Database, start: str, end: str, conf: dict[str, float], model: dict | None = None,
             release_time: bool = True) -> dict:
    """Replay [start, end] one week at a time, as the weekly job runs.

    Each week the automation knows only the mpstruc entries released before
    that week, sees that week's releases, and is scored against what mpstruc's
    curators eventually did with them. release_time=True shows each structure
    as it looked in its release week (no late-linked paper, no later PDB
    annotations); False uses today's PDB data, as a re-check weeks later would.
    Curation lag (curators adding an entry weeks after its release) is not
    modelled; shifting the cutoff by 30-90 days changed the results by <=1 point.
    """
    ids_all = released_between(start, end)
    S = _full_summaries(db)
    kb_full = build_kb(db, S)
    released = {k: (v or {}).get("released", "") for k, v in
                {k: summarise(v) for k, v in rcsb.entry_entities(ids_all).items()}.items()}
    props, stats, truth_ids = [], Counter(), set()
    for ws, we in _weeks(start, end):
        ids = [p for p in ids_all if ws <= released.get(p, "") <= we]
        if not ids:
            continue
        kb = build_kb(db, S, cutoff=ws)
        p, st = propose(kb, ids, conf, model, release_time=release_time)
        props += p
        stats.update(st)
        truth_ids |= {pdb for pdb in ids if pdb in kb_full.by_pdb and pdb not in kb.by_pdb}
    proposed = {p["pdb"]: p for p in props}

    rows = []
    for pdb in sorted(truth_ids | proposed.keys()):
        t = kb_full.by_pdb.get(pdb, [])
        p = proposed.get(pdb)
        row = {"pdb": pdb, "in_mpstruc": bool(t), "proposed": bool(p)}
        if p:
            row.update({"sig": p["sig"], "likelihood": p["likelihood"], "group": p.get("doi") or p.get("pmid") or pdb})
        if t and p:
            tk = t[0]
            row.update({
                "true_kind": tk.kind, "pred_kind": p["placement"],
                "true_sg": [x.sg for x in t], "pred_sg": p["subgroup"],
                "sg_top3": [a["subgroup"] for a in p["alternatives"]],
                "true_owner": _owner(kb_full, tk).pdb if tk.kind == "related" else
                              (_parent_pdb(kb_full, tk.parent_key) if tk.kind == "daughter" else None),
                "pred_owner": p["ownerPdb"], "rule": p["rule"], "confidence": p["confidence"],
            })
        rows.append(row)

    both = [r for r in rows if r["in_mpstruc"] and r["proposed"]]
    def acc(rs, ok):
        return {"n": len(rs), "correct": sum(1 for r in rs if ok(r)),
                "pct": round(100 * sum(1 for r in rs if ok(r)) / len(rs), 1) if rs else None}
    sg_ok = lambda r: r["pred_sg"] in r["true_sg"]
    kind_ok = lambda r: r["pred_kind"] == r["true_kind"]
    full_ok = lambda r: kind_ok(r) and sg_ok(r) and (r["true_kind"] == "parent" or r["pred_owner"] == r["true_owner"])
    by_rule = {rule: {
        "kind": acc([r for r in both if r["rule"] == rule], kind_ok),
        "full": acc([r for r in both if r["rule"] == rule], full_ok),
        "subgroup": acc([r for r in both if r["rule"] == rule and r["pred_sg"]], sg_ok),
        "owner": acc([r for r in both if r["rule"] == rule and r["pred_kind"] in ("daughter", "related", "review")
                      and r["true_kind"] in ("daughter", "related")], lambda r: r["pred_owner"] == r["true_owner"]),
    } for rule in sorted({r["rule"] for r in both})}
    confusion = Counter((r["true_kind"], r["pred_kind"]) for r in both)
    n_true = len(truth_ids)
    weeks = max(1.0, ((dt.date.fromisoformat(end) - dt.date.fromisoformat(start)).days + 1) / 7)
    cuts = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}
    if model and model.get("threshold") is not None:
        cuts.add(model["threshold"])
    bands = []
    for t in sorted(cuts):
        shown = [r for r in rows if r["proposed"] and r["likelihood"] >= t]
        hit = sum(1 for r in shown if r["in_mpstruc"])
        bands.append({"min_likelihood": t, "per_week": round(len(shown) / weeks, 1),
                      "recall_pct": round(100 * hit / n_true, 1) if n_true else None,
                      "precision_pct": round(100 * hit / len(shown), 1) if shown else None})
    auto = [r for r in both if r["pred_kind"] != "review"]
    return {
        "window": [start, end], "release_time": release_time, "released": len(ids_all), "stats": dict(stats),
        "weeks": round(weeks, 1),
        "discovery": {"in_mpstruc": n_true, "found": len(both),
                      "recall_pct": round(100 * len(both) / n_true, 1) if n_true else None,
                      "proposed": len(proposed), "proposed_not_in_mpstruc": len(proposed) - len(both)},
        "bands": bands,
        "placement_type": acc(both, kind_ok),
        "placement_auto": acc(auto, kind_ok),
        "placement_full": acc(both, full_ok),
        "placement_full_auto": acc(auto, full_ok),
        "sent_to_review_pct": round(100 * (len(both) - len(auto)) / len(both), 1) if both else None,
        # Over every structure mpstruc added, including those given no subfamily at all.
        "subgroup_top1": acc(both, sg_ok),
        "subgroup_top3": acc(both, lambda r: bool(set(r["sg_top3"]) & set(r["true_sg"]))),
        "subgroup_suggested_pct": round(100 * sum(1 for r in both if r["pred_sg"]) / len(both), 1) if both else None,
        "by_rule": by_rule,
        "confusion": {f"{a}→{b}": n for (a, b), n in sorted(confusion.items())},
        "training": [(r["sig"], int(r["in_mpstruc"]), r["group"]) for r in rows if r["proposed"]],
        "misses": [r for r in rows if r["in_mpstruc"] and not r["proposed"]][:200],
        "errors": [r for r in both if not full_ok(r)][:300],
    }


def confidences(result: dict) -> dict[str, float]:
    """Per-rule confidence = how often that rule's suggestion was fully right
    (type, subfamily and parent/owner) in the tuning period."""
    out = {}
    names = {"same-paper": "related", "same-paper-batch": "related-batch", "same-protein": "daughter",
             "partial-protein": "partial", "new-protein": "parent"}
    for rule, r in result["by_rule"].items():
        if rule in names and r["full"]["n"]:
            out[names[rule]] = round(r["full"]["correct"] / r["full"]["n"], 3)
    return out


def pick_threshold(training: list, n_true: int, min_recall: float = 0.9, folds: int = 5) -> float:
    """Highest likelihood cut that still finds min_recall of what mpstruc added.

    Uses out-of-fold scores with folds grouped by paper, so structures from
    one paper (which share title words and signals) are never split between
    the fitted and the scored part. n_true counts everything mpstruc added in
    the period, flagged or not.
    """
    import zlib
    fold = lambda g: zlib.crc32(str(g).encode()) % folds
    oof = []
    for f in range(folds):
        train = [r for r in training if fold(r[2]) != f]
        hold = [r for r in training if fold(r[2]) == f]
        m = fit_model([x for x, _, _ in train], [y for _, y, _ in train])
        oof += [(score(x, m), y) for x, y, _ in hold]
    best = 0.0
    for t in [i / 100 for i in range(0, 100)]:
        if sum(y for sc, y in oof if sc >= t) >= min_recall * n_true:
            best = t
    return best


def backtest_md(tune: dict, test: dict, later: dict, conf: dict, model: dict, threshold: float) -> str:
    def headline(r):
        b = next(x for x in r["bands"] if abs(x["min_likelihood"] - threshold) < 1e-9)
        return (f"| {'release week' if r['release_time'] else 'once papers are linked'} | {b['recall_pct']}% | {b['per_week']} | "
                f"{b['precision_pct']}% | {r['subgroup_top1']['pct']}% ({r['subgroup_top3']['pct']}%) | "
                f"{r['placement_auto']['pct']}% | {r['placement_full_auto']['pct']}% | {r['sent_to_review_pct']}% |")

    def block(r, title):
        L = [f"## {title}", "",
             f"Structures released {r['window'][0]} to {r['window'][1]}, replayed one week at a time"
             f"{' as each looked in its release week' if r['release_time'] else ' with today’s PDB data'}.", "",
             f"- Released: {r['released']:,}; later added to mpstruc: {r['discovery']['in_mpstruc']:,}; "
             f"flagged at any likelihood: {r['discovery']['found']:,} ({r['discovery']['recall_pct']}%)", "",
             "| Show if likelihood ≥ | Per week | Finds (recall) | Of those shown, in mpstruc |", "|---|---|---|---|"]
        for b in r["bands"]:
            mark = " ← default" if abs(b["min_likelihood"] - threshold) < 1e-9 else ""
            L.append(f"| {b['min_likelihood']:.2f}{mark} | {b['per_week']} | {b['recall_pct']}% | {b['precision_pct']}% |")
        L += ["",
              f"- Placement type right: {r['placement_type']['pct']}% of all {r['placement_type']['n']}; "
              f"{r['placement_auto']['pct']}% of the {r['placement_auto']['n']} the automation decides "
              f"({r['sent_to_review_pct']}% go to curators with no default)",
              f"- Placement fully right (type + subfamily + parent/owner): {r['placement_full']['pct']}% of all; "
              f"{r['placement_full_auto']['pct']}% of those it decides",
              f"- Subfamily right first time: {r['subgroup_top1']['pct']}%; in top 3: {r['subgroup_top3']['pct']}% "
              f"(of all {r['subgroup_top1']['n']}; a subfamily was suggested for {r['subgroup_suggested_pct']}%)",
              "", "| Rule | Structures | Type right | Fully right | Subfamily right | Parent/owner right |", "|---|---|---|---|---|---|"]
        for rule, x in r["by_rule"].items():
            L.append(f"| {rule} | {x['kind']['n']} | {x['kind']['pct']}% | {x['full']['pct']}% | {x['subgroup']['pct']}% | "
                     f"{str(x['owner']['pct']) + '%' if x['owner']['n'] else '—'} |")
        L += ["", "Confusion (true → predicted): " + ", ".join(f"{k}: {v}" for k, v in r["confusion"].items()), ""]
        return L

    L = ["# Automation backtest", "",
         "Each week of a past period is replayed: the automation knows only the mpstruc entries released before that "
         "week and is scored against what mpstruc's curators eventually did. The likelihood model, the default cut-off "
         "and the per-rule confidences are fitted on Jan–Jun 2025 only (cut-off by cross-validation grouped by paper). "
         "The placement rules themselves were designed from the full 2 Oct 2026 snapshot "
         "(reports/placement_analysis.md), which includes the test period (about 4% of it).", "",
         "## Summary: test period (Jul–Dec 2025, never used for fitting), at the default cut-off", "",
         "| Data the automation has | Finds | Per week | Shown that mpstruc added | Subfamily 1st (top 3) | Type right* | Fully right* | No default |",
         "|---|---|---|---|---|---|---|---|",
         headline(test), headline(later), "",
         "\\* of the structures the automation places itself. *Release week* is what curators see first. The weekly job "
         "re-checks the last 26 weeks every Thursday, so as papers are linked to a structure its suggestion moves toward the "
         "*once papers are linked* row. About half of new releases have no paper linked in their first weeks.", ""]
    L += block(test, "Test period, release week")
    L += block(later, "Test period, once papers are linked")
    L += block(tune, "Tuning period, release week (in-sample: the model was fitted on these)")
    L += ["## Likelihood model (fitted on the tuning period)", "", "| Signal | Weight |", "|---|---|"]
    named = {k: v for k, v in model["weights"].items() if not k.startswith("t:")}
    words = sorted(((k, v) for k, v in model["weights"].items() if k.startswith("t:")), key=lambda kv: -kv[1])
    L += [f"| {_reasons({k: 1.0}, model)[0][0]} | {v:+.2f} |" for k, v in sorted(named.items(), key=lambda kv: -kv[1])]
    L += ["", "Title words that raise the likelihood most: " + ", ".join(f"{k[2:]} ({v:+.2f})" for k, v in words[:15]),
          "", "Title words that lower it most: " + ", ".join(f"{k[2:]} ({v:+.2f})" for k, v in words[::-1][:15]),
          "", "## Confidence shown for each rule (fully right, tuning period)", ""]
    L += [f"- {k}: {round(100 * v)}%" for k, v in conf.items()]
    return "\n".join(L) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--since", help="first release date to consider (default: 26 weeks ago)")
    r.add_argument("--until")
    sub.add_parser("backtest")
    args = ap.parse_args()
    db = load(source_xml())
    SITE_DATA.mkdir(exist_ok=True)

    if args.cmd == "backtest":
        tune = evaluate(db, "2025-01-01", "2025-06-30", {})
        model = fit_model([x for x, _, _ in tune["training"]], [y for _, y, _ in tune["training"]])
        threshold = pick_threshold(tune["training"], tune["discovery"]["in_mpstruc"], min_recall=0.9)
        model["threshold"] = threshold
        tune = evaluate(db, "2025-01-01", "2025-06-30", {}, model)
        conf = confidences(tune)
        test = evaluate(db, "2025-07-01", "2025-12-31", conf, model)
        later = evaluate(db, "2025-07-01", "2025-12-31", conf, model, release_time=False)
        settings = {"confidence": conf, "model": model, "threshold": threshold,
                    "fitted_on": tune["window"], "tested_on": test["window"]}
        (ROOT / "model.json").write_text(json.dumps(settings, indent=1))
        for x in (tune, test, later):
            x.pop("training", None)
        (ROOT / "reports" / "backtest.json").write_text(json.dumps({"tune": tune, "test": test, "later": later, **settings}, indent=1))
        (ROOT / "reports" / "backtest.md").write_text(backtest_md(tune, test, later, conf, model, threshold))
        keep = ("window", "release_time", "released", "weeks", "discovery", "bands", "placement_type", "placement_auto",
                "placement_full", "placement_full_auto", "sent_to_review_pct", "subgroup_top1", "subgroup_top3",
                "subgroup_suggested_pct",
                "by_rule", "confusion")
        site = {"test": {k: test[k] for k in keep} | {"errors": test["errors"][:150], "misses": test["misses"][:100]},
                "later": {k: later[k] for k in keep}, **settings}
        (SITE_DATA / "backtest.json").write_text(json.dumps(site))
        print(backtest_md(tune, test, later, conf, model, threshold))
        return

    today = dt.date.today()
    since = args.since or (today - dt.timedelta(weeks=26)).isoformat()
    until = args.until or today.isoformat()
    settings_path = ROOT / "model.json"
    settings = json.loads(settings_path.read_text()) if settings_path.exists() else {}
    conf, model = settings.get("confidence", {}), settings.get("model")
    kb = build_kb(db, _full_summaries(db))
    ids = released_between(since, until)
    props, stats = propose(kb, ids, conf, model)
    props.sort(key=lambda p: (-p["likelihood"], p["pdb"]))
    snap = (db.meta.get("lastDatabaseEditDate") or db.meta.get("lastNewProteinDate") or "")[:10]
    age = (today - dt.date.fromisoformat(snap)).days if snap else None
    out = {"generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
           "window": [since, until], "snapshot": db.meta.get("timeStamp"), "snapshotAgeDays": age, "stats": stats,
           "threshold": settings.get("threshold", 0.0), "confidence": conf, "candidates": props}
    (SITE_DATA / "candidates.json").write_text(json.dumps(out, separators=(",", ":")))
    print(f"{stats}; {sum(p['likelihood'] >= out['threshold'] for p in props)} at or above the cut-off; "
          f"wrote {SITE_DATA / 'candidates.json'}")


if __name__ == "__main__":
    main()
