"""Find newly released membrane-protein structures and propose mpstruc placements.

    python -m mpstruc.automation run [--since YYYY-MM-DD] [--until YYYY-MM-DD]
    python -m mpstruc.automation backtest
    python -m mpstruc.automation apply-rules   # curator_rules.json on the current candidates, no weekly run
    python -m mpstruc.automation siblings      # add same-paper structures from outside the window, no weekly run

`run` writes ../data/candidates.json for the site's Automation tab.
`apply-rules` re-reads curator_rules.json, recounts barrel strands, rescores
and re-places the candidates a rule touches in ../data/candidates.json, and
rewrites ../data/curator_rules.json and reports/rule_conflicts.md.
`siblings` adds, for every paper of a candidate, the paper's structures that
were released before or after the window and are not in mpstruc (siblings.py), placed
and drafted like every other candidate; `run` does the same every week.
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

The curators' own rules (curator_rules.json, rules.py) are applied after these:
a subfamily with a typical number of barrel strands (autotransporters 12,
Omp85-TpsB 16) is checked against the strands counted in the structure
(barrel.py), and the suggestion is moved when the count says so (apply_rules).
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import re
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from . import rcsb, rules as rules_mod, siblings as siblings_mod
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


_PDB_IN_TEXT = re.compile(r"(?<![A-Za-z0-9])([0-9](?=[A-Z0-9]{0,2}[A-Z])[A-Z0-9]{3})(?![A-Za-z0-9])")


def mentioned_in_descriptions(db: Database) -> dict[str, str]:
    """PDB codes that mpstruc names inside an entry's description text (codes
    with at least one letter, so years and numbers are not taken), mapped to
    the record key of the entry whose description names them."""
    out: dict[str, str] = {}
    for e in db.entries():
        for code in _PDB_IN_TEXT.findall(e.description or ""):
            if code != e.pdbCode.strip().upper():
                out.setdefault(code, e.key)
    return out


def in_mpstruc(kb: KB, mentioned: dict[str, str], ids: list[str]) -> list[dict]:
    """Where each of these PDB releases already is in mpstruc: as an entry
    (parent/daughter), a related code, or only named in a description."""
    out = []
    for pdb in ids:
        known = kb.by_pdb.get(pdb)
        if known:
            k = known[0]
            owner = kb.by_key.get(k.owner_key)
            out.append({"pdb": pdb, "as": k.kind, "key": k.owner_key if k.kind == "related" else k.key,
                        "ownerPdb": owner.pdb if owner and k.kind == "related" else None,
                        "subgroup": k.sg, "subgroupName": kb.sg_names.get(k.sg)})
        elif pdb in mentioned:
            owner = kb.by_key.get(mentioned[pdb])
            out.append({"pdb": pdb, "as": "description", "key": mentioned[pdb],
                        "ownerPdb": owner.pdb if owner else None,
                        "subgroup": owner.sg if owner else None,
                        "subgroupName": kb.sg_names.get(owner.sg) if owner else None})
    return out


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
    why: dict[str, list[tuple[float, str]]] = defaultdict(list)
    for f in sorted(_features(s)):  # in a fixed order, so the sums and the features named are the same every run
        dist = kb.feat_sg.get(f)
        if not dist:
            continue
        total = sum(dist.values())
        spec = 1.0 / math.log2(1 + len(dist))
        for sg, n in dist.items():
            w = weights[f.split(":")[0]] * spec * n / total
            votes[sg] += w
            why[sg].append((-w, f))
    tot = sum(votes.values())
    top = sorted(votes.items(), key=lambda kv: (-kv[1], kv[0]))[:3]
    # The four features that weigh most for the subgroup (then by name).
    return [(sg, v / tot, [f for _, f in sorted(why[sg])[:4]]) for sg, v in top] if tot else []


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
    mine = protein_ids(s)
    for c in paper_keys(s):
        members = [p for p, cr in batch_leads.get(c, []) if cr & mine]  # (in lead order, see batch_leads)
        lead = members[0] if len(members) > 1 else None
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


NEW_ID_FROM = "2025-01-01"  # IDs 1xxx-3xxx released from here on are the PDB's reissued ones (after 9ZZZ)
GENERIC_NAMES = {"", "protein", "uncharacterized protein", "putative uncharacterized protein", "hypothetical protein"}


def lead_key(pdb: str, s: dict | None) -> tuple:
    """Sort key for 'the lowest PDB ID' in mpstruc's sense: the first deposited.

    Classic IDs (1xxx-9xxx) were issued in order, so among them the lowest ID is
    the earliest deposit and the string order is kept. Once 9ZZZ was reached
    (2025) the PDB issued IDs again from 1xxx-3xxx (e.g. 24UG, 28PN), which sort
    below every classic ID but are newer than all of them, and are not issued in
    string order (28PN was deposited before 24UG). Those come after every classic
    ID, ordered by deposit date (then release date, then ID)."""
    pdb = pdb.upper()
    s = s or {}
    rel = s.get("released") or ""
    if pdb[:1] in "123" and rel >= NEW_ID_FROM:
        return (1, s.get("deposited") or rel, rel, pdb)
    return (0, "", "", pdb)


def protein_ids(s: dict) -> set:
    """What 'the same protein' means between structures of one paper: a shared membrane
    UniProt accession, or the same protein name in the same species (strain variants
    of a protein carry different accessions, e.g. Vibrio cholerae NFeoB in 9D8B and 9PS9)."""
    names = s.get("uniprot_names") or {}
    def specific(v: str) -> bool:  # not 'Porin', not a fusion's 'Soluble cytochrome b562,Receptor'
        name, species = v.split("|", 1)
        return bool(species) and name not in GENERIC_NAMES and "," not in name and len(name.split()) >= 3
    named = {f"name:{names[u]}" for u in s["core"] if u in names and specific(names[u])}
    return set(s["core"]) | named


def batch_leads(cands: dict[str, dict]) -> dict[str, list[tuple[str, set]]]:
    """New structures grouped by paper: {citation: [(pdb, proteins), ...]}, in lead order.

    A paper's structures that share a protein (protein_ids) are listed under the
    lowest PDB ID among them, in the sense of lead_key (the first deposited): in
    the snapshot the entry a paper's related structures hang under is the lowest ID
    70% of the time, the best resolution 52%.
    """
    by_cit: dict[str, list[tuple[str, set]]] = defaultdict(list)
    for pdb in sorted(cands, key=lambda p: lead_key(p, cands[p])):
        s = cands[pdb]
        for c in paper_keys(s):
            by_cit[c].append((pdb, protein_ids(s)))
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
            release_time: bool = False, keep: set[str] | None = None) -> tuple[list[dict], dict]:
    """Proposals for the given releases. release_time=True shows each structure
    as it looked in its release week (see as_at_release), for the backtest.
    keep: place exactly these structures (e.g. those likely to be membrane
    proteins) instead of those whose membrane signals fire."""
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
        if (pdb in keep) if keep is not None else in_pool(sig):
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


def add_listings(props: list[dict], db: Database | None = None,
                 summaries: dict[str, dict | None] | None = None) -> None:
    """Attach a draft mpstruc entry to each suggestion (see listing.py).

    `mp` holds the fields in mpstruc raw markup, used both by the XML download
    and (being escaped text) by the page. A structure suggested as the
    lead of other new structures lists them in its description, as mpstruc does.

    What mpstruc already holds is passed along so the draft reads like the
    curators' own text: for a daughter (or a related entry, or one to review)
    the name mpstruc uses for the same protein in that family, and for a
    related entry the entry it goes under with its related structures, so
    `mp.relatedLine` names only what differs from them. When the family gives
    no name (a new parent, say), the name mpstruc uses for entries with exactly
    the same membrane protein(s) anywhere in the database is reused instead
    (from PR #4's known_names; `summaries` are those of the mpstruc entries,
    looked up when needed).
    """
    from . import listing
    db = db or load(source_xml())
    known = {e.pdbCode.upper(): e for e in db.entries()}
    sg_parents = {sg.key: sg.parents for g in db.groups for sg in g.subgroups}
    by_owner = defaultdict(list)
    for p in props:
        if p["placement"] == "related" and p.get("ownerPdb"):
            by_owner[p["ownerPdb"]].append(p["pdb"])
    plan = {}
    for p in props:
        owner = known.get((p.get("ownerPdb") or "").upper())
        family, group = [], []
        if owner and p["placement"] in ("daughter", "related", "review"):
            parent = owner if owner.kind == "parent" else known.get((owner.masterProteinPdbCode or "").upper())
            family = [parent] + parent.daughters if parent else [owner]
        if p["placement"] == "related" and p.get("ownerPdb"):
            lead = p["ownerPdb"].upper()
            group = [lead] + ([r.pdb.upper() for r in owner.related] if owner else [])
            group += [x for x in sorted(by_owner.get(p["ownerPdb"], [])) if x != p["pdb"]]
        plan[p["pdb"]] = (family, group)
    need = {p["pdb"] for p in props} | {e.pdbCode.upper() for f, _ in plan.values() for e in f} | \
           {x for _, g in plan.values() for x in g}
    det = listing.enrich(rcsb.entry_details(sorted(need)))
    by_core: dict[frozenset, list] | None = None

    def same_core() -> dict[frozenset, list]:
        """mpstruc entries by their exact set of membrane proteins (built once, when first needed)."""
        nonlocal by_core, summaries
        if by_core is None:
            by_core = {}
            try:
                summaries = summaries if summaries is not None else _full_summaries(db)
            except Exception as err:  # the drafts still work without it
                print(f"listing: could not look up mpstruc names by protein ({err})")
                summaries = {}
            for e in db.entries():
                s = summaries.get(e.pdbCode.upper())
                if s and s.get("core"):
                    by_core.setdefault(frozenset(s["core"]), []).append(e)
        return by_core

    by_pdb = {p["pdb"]: p for p in props}
    for p in props:
        core = set(p.get("uniprot") or []) - FUSION_UNIPROT
        lead = by_pdb.get(p.get("ownerPdb") or "") if p["placement"] == "related" else None
        if lead:
            # A related structure of a new lead whose extra membrane protein is in no other
            # structure of the study is drafted as its lead's protein, the way the lead's
            # description lists it: 9WFL (GGCX with the membrane protein PRGP1 bound, each of
            # its six siblings with another partner) is 'GGCX with bound PRGP1', not a
            # GGCX-PRGP1 complex. (A KCNQ2/3 heteromer among several stays a heteromer.)
            lead_core = set(lead.get("uniprot") or []) - FUSION_UNIPROT
            mates = [by_pdb[x] for x in by_owner.get(p["ownerPdb"], []) if x != p["pdb"] and x in by_pdb]
            if lead_core < core and len(mates) >= 2 and \
                    all(set(m.get("uniprot") or []) - FUSION_UNIPROT == lead_core for m in mates):
                core = lead_core
        family, group = plan[p["pdb"]]

        def proteins(e):
            d = det.get(e.pdbCode.upper()) or {}
            return {u for x in d.get("polymer_entities") or []
                    for u in (x.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []}
        # Reuse mpstruc's name only from entries of the same protein: all of its
        # membrane proteins if any entry has them, else any of them (a review
        # case shares just part of it, so its family name may not fit).
        names = [(e.name, e.bibliography.year) for e in family if core and core <= proteins(e)] or \
                [(e.name, e.bibliography.year) for e in family if core & proteins(e)]
        if not names and p["placement"] in ("daughter", "related"):
            names = [(e.name, e.bibliography.year) for e in family]
        if not names and core:
            names = [(e.name, e.bibliography.year) for e in same_core().get(frozenset(core), [])]
        ctx = {"family_names": names, "group": [(x, det.get(x)) for x in group],
               "subgroup_names": [e.name for e in sg_parents.get(p.get("subgroup"), [])]}
        mp = listing.build(det.get(p["pdb"]), core, [(r, det.get(r)) for r in sorted(by_owner.get(p["pdb"], []))], ctx)
        if not mp:
            continue
        # Name, species, description, authors and title are built from escaped
        # text plus <br />, <sub>, <sup> and entities only, so the page shows them as HTML.
        p["mp"] = mp

# ---------------------------------------------------------------- curators' rules

PLACEMENT_FIELDS = ("placement", "subgroup", "subgroupName", "group", "ownerPdb", "ownerKey", "confidence", "rule",
                    "evidence", "alternatives")
RULE_FIELDS = ("barrelStrands", "ruleCheck")
MIN_BARREL = 8  # strands (membrane.CORE_MIN_STRANDS): fewer is not a barrel that crosses the membrane
REVIEW_BELOW = 0.5  # membrane score under which family votes alone do not place a structure in a ruled subfamily


def strand_rules(kb: KB) -> dict[str, dict]:
    """{subgroup key: rule} for the curators' rules that give a number of barrel strands."""
    return rules_mod.strand_rules(rules_mod.resolve(rules_mod.load_rules(), kb.sg_names))


def touched_by_rules(c: dict, ruled: dict[str, dict]) -> bool:
    """A candidate whose suggested subfamily, or one of its top subfamily votes, has a rule."""
    return c.get("subgroup") in ruled or any(a.get("subgroup") in ruled for a in c.get("alternatives") or []) \
        or bool(c.get("ruleCheck"))


def _closest_family(kb: KB, sg: str, iprs: set[str]) -> Known | None:
    """The entry of a subfamily whose protein families are most like these (None if none shares any)."""
    best, best_key = None, (0.0, False, "")
    for k in kb.entries:
        if k.sg != sg or k.kind == "related" or not k.s:
            continue
        other = k.s.get("uipr") or set()
        j = len(iprs & other) / len(iprs | other) if iprs | other else 0.0
        key = (round(j, 4), k.kind == "parent", k.released)
        if j > 0 and key > best_key:
            best, best_key = k, key
    return best


def _short(name: str) -> str:
    return name if len(name) <= 60 else name[:57].rstrip() + "..."


UNIPROT_STRANDS = "UniProt strands"  # "method" of a barrel known only from UniProt's transmembrane strands, all resolved


def apply_rules(props: list[dict], kb: KB, ruled: dict[str, dict], barrels: dict[str, dict[str, dict]],
                om: dict[str, set[str]], fams: dict[str, dict | None], conf: dict[str, float],
                conflicts: list[dict] | None = None) -> list[dict]:
    """Check suggestions against the curators' barrel-strand rules, in place.

    props: candidates as propose() returns them, with "membrane" (and "tmEvidence") set.
    ruled: strand_rules. barrels: membrane.counted_barrels output for the touched
    candidates ({pdb: {accession: barrel}}; {} for a structure counted without a barrel;
    a structure that could not be counted is absent); om: {pdb: its outer-membrane
    proteins}; fams: rcsb.entry_families items; conflicts: rule_conflicts output (official
    entries whose own barrel contradicts the rule of the subfamily they are listed in).

    For a candidate whose suggested subfamily or top subfamily votes include a ruled one:
      - "barrelStrands": {"chain", "count", "method", "closed", "uniprot"[, "chains",
        "perChain", "upTo"]} when a membrane barrel was counted in it (or, when none
        could be, UniProt's transmembrane strands of an outer-membrane protein are all
        resolved in it: "method": UNIPROT_STRANDS);
      - the barrel that is checked is the best one known for the paper: a closed barrel
        in the structure itself, else one that decides in its own open barrel, else a
        closed barrel of the same protein in a structure of the same paper ("from"), else
        UniProt's strands. The count of an open barrel is a lower bound (rules.verdict);
      - if the count contradicts the suggested subfamily's rule and fits another ruled
        subfamily, the suggestion moves there: to the entries of the same protein there if
        any, else with the closest entry of the same protein family as closest match (a new
        parent stays a new parent). A structure related to an official entry by its paper
        stays with that entry, and a barrel made by several chains moves nothing; the
        contradiction is only noted;
      - if its own barrel does not decide (none counted, or open with too few strands)
        and the official entry it would join is itself in conflicts, with a barrel of the
        same protein that fits another ruled subfamily, it moves there on that basis;
      - a new protein placed by family votes alone, counted without a barrel and with a
        membrane score under REVIEW_BELOW, goes to review (without a subfamily if none of
        its proteins is an outer-membrane protein at all). A structure that
        could not be counted (absent from barrels) is left as it is ("unchecked");
      - "ruleCheck": {"status": "fits" | "moved" | "differs" | "unclear" | "no-barrel" |
        "unchecked", "text": the same in a sentence, "rule": the ruled subfamily it was
        checked against[, "from", "closest": "family" when the owner is only the closest
        entry of the same protein family, "lines": how many evidence lines the rule added]}.
    Structures of the same paper follow their lead (same-paper-batch), as in propose().
    Returns [{"pdb", "before", "after"}] for the suggestions that changed."""
    from .barrel import METHOD, strands_text
    by_pdb = {c["pdb"]: c for c in props}
    odd = {x["pdb"]: x for x in conflicts or []}
    changes = []
    followers: dict[str, list[dict]] = defaultdict(list)
    for c in props:
        if c.get("rule") == "same-paper-batch" and c.get("ownerPdb"):
            followers[c["ownerPdb"]].append(c)

    def where(c: dict) -> str:
        return f"{c['placement']} in {c.get('subgroupName') or 'no subfamily'}" + (f" ({c['ownerPdb']})" if c.get("ownerPdb") else "")

    def others(sg: str) -> list[dict]:
        return [r for k, r in ruled.items() if k != sg]

    def how(b: dict, sg: str) -> str:
        return rules_mod.verdict(ruled[sg], b["count"], b.get("closed", True), others(sg), b.get("upTo"))

    def counted(pdb: str, sg: str | None) -> list[dict]:
        """The membrane barrels counted in a structure. One in a protein UniProt does not
        place in the outer membrane ("vote") only stands if it is a closed single-chain
        barrel with the number of strands of the suggested subfamily."""
        out = []
        for acc, b in sorted((barrels.get(pdb) or {}).items()):
            if b.get("count", 0) < MIN_BARREL:
                continue
            if b.get("vote") and not (sg in ruled and b.get("closed", True) and not b.get("chains")
                                      and rules_mod.fits(ruled[sg], b["count"])):
                continue
            out.append(dict(b, uniprot=acc))
        return out

    def pick(pdb: str, sg: str | None, accs: set[str] | None = None) -> dict | None:
        """The barrel to check: one that fits the suggested subfamily's rule if any, else the largest."""
        found = [b for b in counted(pdb, sg) if accs is None or b["uniprot"] in accs]
        if not found:
            return None
        fit = [b for b in found if sg in ruled and how(b, sg) == "fits"]
        return max(fit or found, key=lambda b: (b["count"], b["closed"]))

    def uniprot_strands(c: dict) -> dict | None:
        """UniProt's transmembrane strands of an outer-membrane protein, all resolved in the
        structure (from its per-protein table): the barrel, when none could be counted."""
        rows = [r for r in c.get("tmEvidence") or [] if isinstance(r, dict) and r.get("uniprot") in (om.get(c["pdb"]) or ())
                and r.get("status") == "resolved" and isinstance(r.get("segments"), int)
                and r["segments"] >= MIN_BARREL and r.get("resolved") == r["segments"]]
        if not rows:
            return None
        r = max(rows, key=lambda r: r["segments"])
        return {"count": r["segments"], "closed": True, "method": UNIPROT_STRANDS, "uniprot": r["uniprot"]}

    def outcome(b: dict | None, sg: str) -> tuple[str | None, str | None]:
        """("fits" | "unclear" | "differs" | None, the ruled subfamily it fits instead)."""
        if not b:
            return None, None
        v = how(b, sg)
        if v != "differs":
            return v, None
        to = [k for k in ruled if k != sg and how(b, k) == "fits"]
        return "differs", to[0] if to else None

    def evidence_for(c: dict, sg: str) -> dict | None:
        """The best barrel known for a lead and the structures of its paper."""
        own = pick(c["pdb"], sg)
        decides = lambda b: outcome(b, sg)[0] == "fits" or bool(outcome(b, sg)[1])
        if own and (own.get("closed", True) or decides(own)):
            return own
        accs = {own["uniprot"]} if own else set(c.get("uniprot") or [])
        mates = [dict(b, **{"from": f["pdb"]}) for f in followers.get(c["pdb"], [])
                 for b in [pick(f["pdb"], sg, accs)] if b and b.get("closed", True)]
        if mates:
            return max(mates, key=lambda b: (how(b, sg) == "fits", b["count"]))
        return uniprot_strands(c) or own

    def note_barrel(c: dict, b: dict | None) -> None:
        if b and not b.get("from"):
            c["barrelStrands"] = {k: b[k] for k in ("chain", "count", "method", "closed", "uniprot", "chains", "perChain", "upTo")
                                  if b.get(k) is not None}
        elif c["pdb"] in barrels:
            c.pop("barrelStrands", None)

    def chain_text(b: dict) -> str:
        return (f"{b['chains']} chains together" if b.get("chains") else f"chain {b['chain']}") + \
            ("" if b.get("closed", True) else ", seam open")

    def counted_text(b: dict) -> str:
        """'16 strands counted in the barrel (chain A)' and its variants."""
        if b.get("method") == UNIPROT_STRANDS:
            return (f"all {b['count']} membrane-crossing strands UniProt gives for {b['uniprot']} are resolved "
                    f"(no barrel could be counted from the coordinates)")
        return (f"{strands_text(b)} strands counted in the barrel ({chain_text(b)}"
                + (f" of {b['from']}, same paper and protein" if b.get("from") else "") + ")")

    def annotate(c: dict, status: str, rule: str, text: str, lines: list[str] | None = None, **more) -> None:
        text = text[0].upper() + text[1:]
        lines = [text] + (lines or [])
        c["evidence"] = list(c["evidence"]) + lines
        c["ruleCheck"] = {"status": status, "rule": rule, "text": text, "lines": len(lines), **more}

    def move(c: dict, cur: str, to: str, text: str, acc: str | None) -> None:
        rule, new = ruled[cur], ruled[to]
        ev = []
        if c.get("ownerPdb") and c["ownerPdb"] not in text:
            ev.append(f"mpstruc lists the same protein under {rule['name']} ({c['ownerPdb']}); "
                      f"by the rule that entry belongs with {new['name']} too")
        core = set(c.get("uniprot") or []) - FUSION_UNIPROT
        family = False  # the closest match is an entry of the same family, not of the same protein
        same = [k for u in core for k in kb.by_uniprot.get(u, []) if k.sg == to]
        exact = [k for k in same if k.s and k.s["core"] == core]
        if exact:
            votes_p = Counter(k.parent_key for k in exact)
            parent_key = sorted(votes_p, key=lambda key: (-votes_p[key], key))[0]
            c.update(placement="daughter", ownerPdb=_parent_pdb(kb, parent_key), ownerKey=parent_key,
                     confidence=round(conf.get("daughter", 0.9), 2))
            ev.append(f"Same protein(s) as {votes_p[parent_key]} structure(s) under {c['ownerPdb']} ({parent_key}) there")
        elif same:
            parent_key = Counter(k.parent_key for k in same).most_common(1)[0][0]
            c.update(placement="review", ownerPdb=_parent_pdb(kb, parent_key), ownerKey=parent_key,
                     confidence=round(conf.get("partial", 0.6), 2))
            ev.append(f"Shares a protein with {c['ownerPdb']} ({parent_key}) there (complex or fusion differs)")
        else:
            iprs = set(((fams.get(c["pdb"]) or {}).get("ipr") or {}).get(acc) or ())
            near = _closest_family(kb, to, iprs)
            if c["placement"] == "parent":
                c.update(ownerPdb=None, ownerKey=None)
                if near:
                    ev.append(f"Closest entry there: {near.pdb} (same protein family)")
            else:  # it was a daughter of, or shared a protein with, an entry of the other subfamily
                family = bool(near)
                c.update(placement="review", ownerPdb=_parent_pdb(kb, near.parent_key) if near else None,
                         ownerKey=near.parent_key if near else None, confidence=0.0)
                if near:
                    ev.append(f"Closest entry there: {near.pdb} (same protein family, another species or complex): "
                              f"new parent, or daughter of it")
        c.update(subgroup=to, subgroupName=kb.sg_names.get(to), group=kb.sg_group.get(to), rule="barrel-rule")
        annotate(c, "moved", to, text, ev, **{"from": cur}, **({"closest": "family"} if family else {}))

    for c in props:
        c.pop("ruleCheck", None)
    leads = [c for c in props if c.get("rule") != "same-paper-batch"]
    for c in leads:
        if not touched_by_rules(c, ruled):
            continue
        before = where(c)
        cur = c.get("subgroup")
        own = pick(c["pdb"], cur)
        b = evidence_for(c, cur) if cur in ruled else own
        note_barrel(c, own or (b if b and b.get("method") == UNIPROT_STRANDS else None))
        if cur not in ruled:
            continue
        rule = ruled[cur]
        status, to = outcome(b, cur)
        several = bool(b and b.get("chains"))
        # The official entry it would join is itself listed against the rule (TamA 4C00).
        o = odd.get(c.get("ownerPdb") or "")
        by_owner = None
        if o and o.get("subgroup") == cur and o.get("fits") in ruled and c["placement"] != "related" \
                and o.get("uniprot") in (c.get("uniprot") or []) and status != "fits" and not (status == "differs" and to) \
                and (b is None or rules_mod.possible(ruled[o["fits"]], b["count"], b.get("closed", True), b.get("upTo"))):
            by_owner = o["fits"]
        if status == "fits":
            annotate(c, "fits", cur, f"{counted_text(b)}: fits the curators' rule for {rule['name']} ({rules_mod.strands_phrase(rule)})")
        elif status == "differs" and to and c["placement"] != "related" and not several:
            new = ruled[to]
            move(c, cur, to, f"{counted_text(b)}: {new['name']} are {rules_mod.strands_phrase(new)}, {rule['name']} "
                             f"{rules_mod.strands_phrase(rule)} (curators' rule), so it is suggested under {new['name']}", b["uniprot"])
        elif by_owner:
            new = ruled[by_owner]
            own_part = "" if b is None else f"; in this structure, {counted_text(b)}"
            move(c, cur, by_owner,
                 f"mpstruc lists the same protein under {rule['name']} ({o['pdb']}), but "
                 f"{o['count']} strands are counted in the barrel of {o['pdb']} itself: {new['name']} are "
                 f"{rules_mod.strands_phrase(new)}, {rule['name']} {rules_mod.strands_phrase(rule)} (curators' rule), so it is "
                 f"suggested under {new['name']}{own_part}", o["uniprot"])
        elif status == "differs":
            text = f"{counted_text(b)}, but {rule['name']} are {rules_mod.strands_phrase(rule)} (curators' rule)"
            if to and c["placement"] == "related":
                text += (f"; it fits {ruled[to]['name']} ({rules_mod.strands_phrase(ruled[to])}). It stays with "
                         f"{c.get('ownerPdb')}, the official entry of the same paper")
            elif to:
                text += (f"; the number fits {ruled[to]['name']} ({rules_mod.strands_phrase(ruled[to])}), but those barrels are "
                         f"made by one chain, so it is not moved")
            else:
                may = [r for r in others(cur) if rules_mod.possible(r, b["count"], b.get("closed", True), b.get("upTo"))]
                if may and not b.get("closed", True):
                    text += (f"; the barrel is open, so strands may be missing from the count: it may be one of "
                             f"{may[0]['name']} ({rules_mod.strands_phrase(may[0])})")
            annotate(c, "differs", cur, text)
        elif status == "unclear":
            if b.get("closed", True):  # two rings fit the same strands
                text = (f"{counted_text(b)}: only one of the two numbers fits the curators' rule for {rule['name']} "
                        f"({rules_mod.strands_phrase(rule)}), so the count does not decide and the suggestion stays")
            elif rules_mod.fits(rule, b["count"]):
                more = [r for r in others(cur) if r["barrelStrands"] > rule["barrelStrands"]
                        and rules_mod.possible(r, b["count"], False)]
                text = (f"{counted_text(b)}: that would fit {rule['name']} ({rules_mod.strands_phrase(rule)}, curators' rule), but "
                        f"the barrel is open, so strands may be missing from the count, and {more[0]['name']} are "
                        f"{rules_mod.strands_phrase(more[0])}: the count does not decide and the suggestion stays")
            else:
                text = (f"{counted_text(b)}, fewer than {rule['name']} typically have ({rule['barrelStrands']}, curators' rule); "
                        f"the barrel is open, so strands may be missing from the count and the suggestion stays")
            annotate(c, "unclear", cur, text)
        elif c["pdb"] not in barrels:
            # Not "no barrel": its chains could not be looked at this time (an outer-membrane
            # protein, or any protein of a small structure suggested here: rule_inputs).
            c["ruleCheck"] = {"status": "unchecked", "rule": cur, "lines": 0,
                              "text": f"The barrel could not be counted in this run (the coordinates could not be fetched), so the "
                                      f"curators' rule for {rule['name']} ({rules_mod.strands_phrase(rule)}) was not checked"}
        elif c.get("rule") == "new-protein" and (c.get("membrane") or 0) < REVIEW_BELOW:
            has_om = bool(om.get(c["pdb"]))
            text = (f"No membrane beta-barrel is resolved in this structure (membrane protein {round(100 * (c.get('membrane') or 0))}% "
                    f"likely), and only its protein families point to {rule['name']}, a barrel subfamily "
                    f"({rules_mod.strands_phrase(rule)}, curators' rule): "
                    + ("only a part outside the membrane of an outer-membrane protein is in it, so a curator decides whether to list it"
                       if has_om else "none of its proteins is an outer-membrane barrel protein, so it is not placed there"))
            c.update(placement="review", confidence=0.0, rule="barrel-rule-no-barrel", ownerPdb=None, ownerKey=None)
            if not has_om:
                c.update(subgroup=None, subgroupName=None, group=None)
            annotate(c, "no-barrel", cur, text)
        if where(c) != before:
            changes.append({"pdb": c["pdb"], "before": before, "after": where(c)})
    for c in props:
        if c.get("rule") != "same-paper-batch":
            continue
        lead = by_pdb.get(c.get("ownerPdb") or "")
        if not lead or not (lead.get("ruleCheck") or touched_by_rules(c, ruled)):
            continue
        before = where(c)
        sg = lead.get("subgroup")
        b = pick(c["pdb"], sg) or uniprot_strands(c)
        note_barrel(c, b)
        check = lead.get("ruleCheck")
        if check and check["status"] in ("moved", "no-barrel"):
            c.update(subgroup=lead.get("subgroup"), subgroupName=lead.get("subgroupName"), group=lead.get("group"))
            text = (f"Follows {lead['pdb']} (same paper), which the curators' barrel rule "
                    + (f"moved to {lead['subgroupName']}" if check["status"] == "moved" else "sent to review")
                    + (f"; in this structure, {counted_text(b)}" if b else ""))
            annotate(c, check["status"], check["rule"], text, **({"from": check["from"]} if check.get("from") else {}))
        elif b and c.get("subgroup") in ruled:
            rule = ruled[c["subgroup"]]
            v = how(b, c["subgroup"])
            text = (f"{counted_text(b)}: fits the curators' rule for {rule['name']} ({rules_mod.strands_phrase(rule)})" if v == "fits"
                    else f"{counted_text(b)}: the barrel is open, so strands may be missing from the count, which does not "
                         f"decide between the curators' rules ({rule['name']}: {rules_mod.strands_phrase(rule)})" if v == "unclear"
                    else f"{counted_text(b)}, but {rule['name']} are {rules_mod.strands_phrase(rule)} (curators' rule)")
            annotate(c, v, c["subgroup"], text)
        if where(c) != before:
            changes.append({"pdb": c["pdb"], "before": before, "after": where(c)})
    return changes


def rule_inputs(props: list[dict], ruled: dict[str, dict]) -> tuple[dict, dict, dict]:
    """(barrels, outer-membrane proteins, families) apply_rules needs, for the candidates
    a rule touches (and the structures of the same paper that follow them). A structure
    whose barrel could not be counted is absent from barrels (if the coordinate server
    cannot be reached, all but those with no chain to look at: membrane.nothing_to_count)."""
    from . import membrane
    lead_touched = {c["pdb"] for c in props if touched_by_rules(c, ruled)}
    ids = sorted(lead_touched | {c["pdb"] for c in props
                                 if c.get("rule") == "same-paper-batch" and c.get("ownerPdb") in lead_touched})
    if not ids:
        return {}, {}, {}
    ext, fams, ents = rcsb.entry_extras(ids), rcsb.entry_families(ids), rcsb.entry_entities(ids)
    om = {p: membrane.outer_membrane_proteins((ext.get(p) or {}).get("kw"), (fams.get(p) or {}).get("ipr")) for p in ids}
    # Suggested under a ruled subfamily: the vote is the membrane evidence, so every chain
    # is looked at when UniProt names no outer-membrane protein in it (invasin).
    every = {c["pdb"] for c in props if c.get("subgroup") in ruled and c["pdb"] in ids}
    try:
        barrels = membrane.counted_barrels(ids, ext, fams, None, ents, every=every)
    except Exception as err:  # without counts the rules only note what they could not check
        print(f"rules: barrels could not be counted ({err})")
        barrels = {p: {} for p in ids if membrane.nothing_to_count(p, ext, fams, ents, every)}
    return barrels, om, fams


def rule_conflicts(db: Database, kb: KB, ruled: dict[str, dict]) -> list[dict]:
    """Official mpstruc entries (parents and daughters) in a ruled subfamily whose counted
    barrel contradicts the rule: information for the curators, mpstruc itself is not changed.

    An entry with several barrels (BAM with a substrate) is fine if one fits, and an open
    barrel with too few strands says nothing (rules.verdict). Entries with no membrane barrel
    in the structure (POTRA, passenger or head domains, BamB-E) are not listed, and neither
    are entries whose protein UniProt does not place in the outer membrane (invasin): their
    chains are not looked at here (membrane.counted_barrels)."""
    from . import membrane
    entries = [(sg.key, e) for g in db.groups for sg in g.subgroups if sg.key in ruled
               for p in sg.parents for e in (p, *p.daughters)]
    ids = sorted({e.pdbCode.upper() for _, e in entries})
    ext, fams, ents = rcsb.entry_extras(ids), rcsb.entry_families(ids), rcsb.entry_entities(ids)
    barrels = membrane.counted_barrels(ids, ext, fams, None, ents)
    out = []
    for sg, e in entries:
        pdb = e.pdbCode.upper()
        found = [dict(b, uniprot=acc) for acc, b in sorted((barrels.get(pdb) or {}).items()) if b.get("count", 0) >= MIN_BARREL]
        rest = lambda k: [r for x, r in ruled.items() if x != k]
        how = lambda b, k: rules_mod.verdict(ruled[k], b["count"], b["closed"], rest(k), b.get("upTo"))
        if not found or any(how(b, sg) != "differs" for b in found):
            continue
        b = max(found, key=lambda b: b["count"])
        other = [k for k in ruled if k != sg and how(b, k) == "fits"]
        out.append({"pdb": pdb, "key": e.key, "name": to_text(e.name) or "", "subgroup": sg, "subgroupName": ruled[sg]["name"],
                    "expected": ruled[sg]["barrelStrands"], "count": b["count"], "chain": b["chain"], "closed": b["closed"],
                    "uniprot": b["uniprot"], "fits": other[0] if other else None,
                    "fitsName": ruled[other[0]]["name"] if other else None})
    return out


def conflicts_md(conflicts: list[dict], ruled: dict[str, dict]) -> str:
    L = ["# Official mpstruc entries that contradict a curators' rule", "",
         "Barrel strands counted from the coordinates (`mpstruc/barrel.py`) against the rules in `curator_rules.json`.",
         "For the curators' information: the automation cannot change official mpstruc. Regenerate with",
         "`python -m mpstruc.automation apply-rules`.", "", "## Rules", ""]
    L += [f"- **{r['name']}** ({k}): {r['note']} (added by {r.get('addedBy', '?')}, {r.get('addedOn', '?')}; "
          f"counts within {r.get('tolerance', 0)} of {r['barrelStrands']} fit)" for k, r in sorted(ruled.items())]
    L += ["", f"## Entries that do not fit ({len(conflicts)})", ""]
    if conflicts:
        L += ["| PDB | Entry | Listed under | Rule | Counted | Fits instead |", "|---|---|---|---|---|---|"]
        for c in conflicts:
            L.append(f"| {c['pdb']} | {c['name']} ({c['key']}) | {c['subgroupName']} | {c['expected']} strands | "
                     f"{'' if c['closed'] else 'at least '}{c['count']} strands (chain {c['chain']}{'' if c['closed'] else ', seam open'}) | "
                     f"{c['fitsName'] or 'no ruled subfamily'} |")
    else:
        L.append("None.")
    return "\n".join(L) + "\n"


def _conflicts(db: Database, kb: KB, ruled: dict[str, dict]) -> list[dict] | None:
    """rule_conflicts, or None when the official entries could not be counted this time."""
    try:
        return rule_conflicts(db, kb, ruled) if ruled else []
    except Exception as err:
        print(f"rules: the official entries could not be checked ({err})")
        return None


def write_rule_files(db: Database, kb: KB, generated: str, conflicts: list[dict] | None = None) -> dict:
    """../data/curator_rules.json for the site (the rules, the official entries that
    contradict one, how well the counting did) and reports/rule_conflicts.{md,json}.
    conflicts: rule_conflicts output if already computed."""
    resolved = rules_mod.resolve(rules_mod.load_rules(), kb.sg_names)
    ruled = rules_mod.strand_rules(resolved)
    if conflicts is None:
        conflicts = rule_conflicts(db, kb, ruled) if ruled else []
    val_path = ROOT / "reports" / "barrel_validation.json"
    validation = None
    if val_path.exists():
        v = json.loads(val_path.read_text())
        validation = {"n": v["n"], "correct": v["correct"], "withinOne": v["within_one"]}
        if v.get("open_n"):
            validation.update(openN=v["open_n"], openCorrect=v["open_correct"])
    site = rules_mod.site_data(resolved, conflicts, generated, validation)
    SITE_DATA.mkdir(exist_ok=True)
    (SITE_DATA / "curator_rules.json").write_text(json.dumps(site, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    (ROOT / "reports" / "rule_conflicts.md").write_text(conflicts_md(conflicts, ruled), encoding="utf-8")
    (ROOT / "reports" / "rule_conflicts.json").write_text(json.dumps(conflicts, indent=1, ensure_ascii=False), encoding="utf-8")
    return site


def reapply_rules(path: Path | None = None, db: Database | None = None) -> dict:
    """Apply curator_rules.json to data/candidates.json without a weekly run.

    1. membrane.backfill_table: barrels are recounted, and the membrane score, reasons and
       per-protein table of the structures a counted barrel changes are replaced.
    2. The candidates a rule touches are placed again from scratch (same classification as
       the weekly run, from the caches) and checked against the rules (apply_rules); their
       draft listing is rebuilt when the placement changed.
    Everything else in the file stays as it is; it is written as the weekly run writes it.
    Returns {"scores": [...], "placements": [...], "reclassified": [...]} (before -> after)."""
    from . import membrane
    path = path or SITE_DATA / "candidates.json"
    db = db or load(source_xml())
    scores = membrane.backfill_table(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    cands = data.get("candidates") or []
    settings_path = ROOT / "model.json"
    settings = json.loads(settings_path.read_text()) if settings_path.exists() else {}
    conf, model = settings.get("confidence", {}), settings.get("model")
    kb = build_kb(db, _full_summaries(db))
    ruled = strand_rules(kb)
    touched = {c["pdb"] for c in cands if touched_by_rules(c, ruled)}
    touched |= {c["pdb"] for c in cands if c.get("rule") == "same-paper-batch" and c.get("ownerPdb") in touched}
    out = {"scores": scores["score_changed"], "placements": [], "reclassified": [], "touched": len(touched)}
    conflicts = None
    if touched:
        ids = [c["pdb"] for c in cands]
        fresh, _ = propose(kb, ids, conf, model, keep=set(ids))
        fresh = {p["pdb"]: p for p in fresh}
        work = []
        for c in cands:
            f = fresh.get(c["pdb"])
            if c["pdb"] not in touched or not f:
                continue
            w = {**c, **{k: f[k] for k in PLACEMENT_FIELDS}}
            w.pop("ruleCheck", None)
            # The same classification as stored: keep the stored wording, without the lines a
            # rule added last time ("lines"). A candidate a rule moved gets the fresh wording,
            # which is the same every run (subgroup_scores names the features in a fixed order).
            same = all(c.get(k) == f[k] for k in ("placement", "subgroup", "ownerPdb", "ownerKey", "rule", "confidence"))
            votes = lambda alts: [(a.get("subgroup"), a.get("score")) for a in alts or []]
            added = (c.get("ruleCheck") or {}).get("lines", None if c.get("ruleCheck") else 0)
            if same and isinstance(added, int) and 0 <= added <= len(c.get("evidence") or []):
                w["evidence"] = list(c.get("evidence") or [])[:len(c.get("evidence") or []) - added]
            if votes(c.get("alternatives")) == votes(f["alternatives"]):
                w["alternatives"] = c["alternatives"]
            work.append(w)
        # A lead that is itself not touched keeps its stored suggestion; its followers need it to follow.
        known = {w["pdb"] for w in work}
        context = work + [c for c in cands if c["pdb"] not in known]
        barrels, om, fams = rule_inputs(work, ruled)
        conflicts = _conflicts(db, kb, ruled)
        apply_rules(context, kb, ruled, barrels, om, fams, conf, conflicts)
        by_pdb = {c["pdb"]: c for c in cands}
        show = lambda c: f"{c['placement']} in {c.get('subgroupName') or 'no subfamily'}" + (f" ({c['ownerPdb']})" if c.get("ownerPdb") else "")
        moved = []
        for w in work:
            c = by_pdb[w["pdb"]]
            f = fresh[w["pdb"]]
            if any(c.get(k) != f[k] for k in ("placement", "subgroup", "ownerPdb")) and not c.get("ruleCheck"):
                out["reclassified"].append({"pdb": c["pdb"], "before": show(c), "after": show(f)})
            if any(c.get(k) != w.get(k) for k in ("placement", "subgroup", "ownerPdb")):
                out["placements"].append({"pdb": c["pdb"], "before": show(c), "after": show(w)})
                moved.append(c["pdb"])
            for k in PLACEMENT_FIELDS:
                c[k] = w[k]
            for k in RULE_FIELDS:
                if k in w:
                    c[k] = w[k]
                else:
                    c.pop(k, None)
        if moved:  # the draft entry names what it goes under: rebuild it for those that moved
            copy = json.loads(json.dumps(cands))
            add_listings(copy, db)
            new_mp = {c["pdb"]: c.get("mp") for c in copy}
            for pdb in moved:
                if new_mp.get(pdb):
                    by_pdb[pdb]["mp"] = new_mp[pdb]
    data["rulesApplied"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes")
    path.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    write_rule_files(db, kb, data["rulesApplied"], conflicts)
    return out


def relist(path: Path | None = None, db: Database | None = None) -> dict:
    """Re-draft the mpstruc listing (`mp`: name, description, relatedLine,
    ligandNotes, ...) of every candidate in data/candidates.json with the
    current listing.py, from the caches (the PDB is asked only for what they
    lack). Nothing is re-scored or re-placed: every other field of the file,
    and its compact JSON format, stay as they are. A candidate whose PDB
    details cannot be had keeps its old draft. Returns what changed."""
    path = path or SITE_DATA / "candidates.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    cands = doc.get("candidates") or []
    # add_listings only reads pdb, placement, ownerPdb and uniprot; drafts go to copies.
    fresh = [{k: p.get(k) for k in ("pdb", "placement", "ownerPdb", "uniprot")} for p in cands]
    add_listings(fresh, db)
    out = {"candidates": len(cands), "redrafted": 0, "names": 0, "descriptions": 0, "relatedLines": 0, "changed": []}
    for p, f in zip(cands, fresh):
        new = f.get("mp")
        if not new:
            continue
        old = p.get("mp") or {}
        merged = {**old, **new}
        if "relatedLine" not in new:
            merged.pop("relatedLine", None)
        out["redrafted"] += 1
        diff = [k for k in ("name", "description", "relatedLine") if old.get(k) != merged.get(k)]
        out["names"] += "name" in diff
        out["descriptions"] += "description" in diff
        out["relatedLines"] += "relatedLine" in diff
        if diff:
            out["changed"].append({"pdb": p["pdb"], **{k: [old.get(k), merged.get(k)] for k in diff}})
        p["mp"] = merged
    path.write_text(json.dumps(doc, separators=(",", ":")), encoding="utf-8")
    return out

# ---------------------------------------------------------------- same-paper siblings

DROPPED_FIELDS = ("sig", "taxa", "likelihood", "reasons")  # of propose(): not written to the site's file
MEMBRANE_FLOOR = 0.1


def _cores(ids: list[str]) -> dict[str, set]:
    """{pdb: its membrane protein(s)} as summarise() gives them (from the caches)."""
    ents, extras, fams = rcsb.entry_entities(ids), rcsb.entry_extras(ids), rcsb.entry_families(ids)
    out = {}
    for p in ids:
        s = summarise(ents.get(p), extras.get(p), fams.get(p))
        if s:
            out[p] = set(s["core"])
    return out


def window_candidates(ids: list[str], fresh_cits: dict[str, dict | None] | None = None) -> list[dict]:
    """[{"pdb", "doi", "pmid", "awaiting"}] for the structures the weekly run keeps, before
    they are placed: their paper from the PDB (the citation just fetched for those tracked
    while they wait, else the cached summary), and whether it is out (awaiting.py)."""
    from . import awaiting
    summ = rcsb.entry_summaries(ids)
    out = []
    for p in ids:
        cit = (fresh_cits or {})[p] if p in (fresh_cits or {}) else (summ.get(p) or {}).get("rcsb_primary_citation")
        cit = cit or {}
        out.append({"pdb": p, "doi": cit.get("pdbx_database_id_DOI"), "pmid": cit.get("pdbx_database_id_PubMed"),
                    "awaiting": bool(awaiting.citation_status(cit or None))})
    return out


def find_siblings(cands: list[dict], window: list[str], have: set[str], kb: KB, mentioned: dict[str, str],
                  decided: set[str], db: Database | None = None, floor: float = MEMBRANE_FLOOR,
                  discover=None, score=None, cores=None,
                  fallback: dict[str, list[str]] | None = None) -> tuple[dict[str, dict], dict[str, dict], dict]:
    """A paper's other structures released before or after the window (siblings.py), scored.

    cands: the window's candidates ({"pdb", "doi", "pmid"[, "awaiting", "uniprot"]});
    have: every structure the run already has;
    fallback: {paper: [pdb, ...]} from the previous run, for a lookup that fails. Returns (siblings kept: {pdb: {"paper",
    "with", "released", "kept": "score" | "same protein as X"}}, their
    membrane.score_structures records, counts). discover, score and cores stand in for
    siblings.discover, membrane.score_structures and _cores in the tests."""
    from . import membrane
    listed = set(kb.by_pdb) | set(mentioned)
    sibs, counts = ((discover or siblings_mod.discover)(cands, window, have, listed, decided, fallback=fallback)
                    if fallback else (discover or siblings_mod.discover)(cands, window, have, listed, decided))
    if not sibs:
        return {}, {}, {**counts, "kept": 0, "belowFloor": 0}
    ids = sorted(sibs)
    mem = score(ids) if score else membrane.score_structures(ids, db)
    kept_cands = {c["pdb"] for c in cands}
    mates = sorted({m for s in sibs.values() for m in s["with"]} & kept_cands)
    core = (cores or _cores)(ids + mates)
    why = siblings_mod.keep(sibs, {p: v["p"] for p, v in mem.items()}, floor, core, kept_cands)
    out = {p: {**sibs[p], "kept": why[p]} for p in ids if p in why}
    counts = {**counts, "kept": len(out), "belowFloor": sum(v["p"] < floor for p, v in mem.items() if p not in why),
              "keptBySameProtein": sum(w != "score" for w in why.values()),
              "papersWithSiblings": len({s["paper"] for s in out.values()})}
    return out, {p: mem[p] for p in out}, counts


def _with_membrane(p: dict, m: dict) -> dict:
    """A propose() suggestion as the site's file holds it."""
    p = {k: v for k, v in p.items() if k not in DROPPED_FIELDS}
    p["membrane"], p["membraneWhy"] = m["p"], m["why"]
    if m.get("tmTable"):
        p["tmEvidence"] = m["tmTable"]
    return p


BATCH_TEXT = ", released in the same batch"  # classify(), rule same-paper-batch


def mark_siblings(props: list[dict], sibs: dict[str, dict], window: list[str]) -> None:
    """Set "sibling" on the suggestions that are siblings, and say what is true of them
    where classify() wrote 'released in the same batch' of a structure and its lead."""
    for p in props:
        own, lead = p["pdb"] in sibs, p.get("ownerPdb") in sibs
        if own:
            p["sibling"] = siblings_mod.note(sibs[p["pdb"]], window)
        if p.get("rule") == "same-paper-batch" and (own or lead) and p.get("evidence"):
            first = str(p["evidence"][0])
            if first.endswith(BATCH_TEXT):
                where = lambda pdb: {"before": "before the window", "after": "after the window",
                                     "inside": "in the window"}[siblings_mod.when(sibs[pdb].get("released"), window)]
                if lead and not own:
                    what = f" (released {sibs[p['ownerPdb']]['released']}, {where(p['ownerPdb'])})"
                elif lead and where(p["pdb"]) == where(p["ownerPdb"]):
                    what = f"; both were released {where(p['pdb'])}"
                elif lead:
                    what = f"; this one was released {where(p['pdb'])}, {p['ownerPdb']} {where(p['ownerPdb'])}"
                else:
                    what = f"; this one was released {where(p['pdb'])}"
                p["evidence"] = [first[:-len(BATCH_TEXT)] + what] + list(p["evidence"][1:])


def add_siblings(path: Path | None = None, db: Database | None = None, find=None) -> dict:
    """Add a paper's structures from outside the window to data/candidates.json, in place.

    For every paper of a candidate whose paper is out, the paper's other PDB entries are
    looked up (siblings.py). Those released before or after the window, not in mpstruc
    and not decided by a curator are scored, placed together with the paper's candidates (same
    code as the weekly run: propose, apply_rules), given their draft listing, and added
    with "sibling": {...}. The paper's candidates already in the file are placed again
    with them (the lead may now be an older structure: batch_leads takes the first deposited),
    keep every field the placement does not own, and have their listing re-drafted (the
    lead's description gains the related lines). Every other candidate, and the compact
    JSON format, stay as they are. Returns {"added": [...], "studies": [{paper, before,
    after}], "counts": {...}}; nothing is written when there is nothing to add."""
    from . import autoreject
    path = path or SITE_DATA / "candidates.json"
    db = db or load(source_xml())
    doc = json.loads(path.read_text(encoding="utf-8"))
    cands = doc.get("candidates") or []
    window = doc["window"]
    settings_path = ROOT / "model.json"
    settings = json.loads(settings_path.read_text()) if settings_path.exists() else {}
    conf, model = settings.get("confidence", {}), settings.get("model")
    kb = build_kb(db, _full_summaries(db))
    mentioned = mentioned_in_descriptions(db)
    curation_path = path.parent / "curation.json"
    curation = json.loads(curation_path.read_text(encoding="utf-8")) if curation_path.exists() else {}
    decided = autoreject.decided_ids(curation)
    # As in the weekly run: a structure a curator decided is not a candidate any more, so
    # it brings in no siblings (it is still in the file until the next run).
    live = [c for c in cands if c["pdb"].upper() not in decided]
    have = {c["pdb"].upper() for c in cands}
    sibs, mem, counts = find_siblings(live, window, have, kb, mentioned, decided, db, discover=find,
                                      fallback=siblings_mod.known_entries(cands))
    out = {"added": [], "studies": [], "counts": counts}
    if not sibs:
        return out
    key_of = lambda c: siblings_mod.paper_key(c.get("doi"), c.get("pmid"))
    touched_papers = {s["paper"] for s in sibs.values()}
    old = {c["pdb"]: c for c in cands if key_of(c) in touched_papers and not c.get("awaiting")}
    ids = sorted(set(old) | set(sibs))
    fresh, _ = propose(kb, ids, conf, model, keep=set(ids))
    fresh = {p["pdb"]: p for p in fresh}
    show = lambda c: f"{c['pdb']} {c['placement']}" + (f" of {c['ownerPdb']}" if c.get("ownerPdb") else "")
    before = defaultdict(list)
    for c in old.values():
        before[key_of(c)].append(show(c))
    work, moved = [], set()
    for pdb in ids:
        f = fresh.get(pdb)
        if not f:
            continue  # the PDB gave nothing to place: an old candidate stays as it is
        if pdb in sibs:
            w = _with_membrane(f, mem[pdb])
        else:
            c = old[pdb]
            w = json.loads(json.dumps(c))
            if any(c.get(k) != f[k] for k in ("placement", "subgroup", "ownerPdb", "ownerKey", "rule")):
                moved.add(pdb)
                w.update({k: f[k] for k in PLACEMENT_FIELDS})
                for k in RULE_FIELDS:
                    w.pop(k, None)
        work.append(w)
    mark_siblings([w for w in work if w["pdb"] in sibs or w["pdb"] in moved], sibs, window)
    # The curators' rules, for the new and the moved ones only (the others keep what the
    # last run or apply-rules noted).
    ruled = strand_rules(kb)
    check = [w for w in work if w["pdb"] in sibs or w["pdb"] in moved]
    if ruled and any(touched_by_rules(w, ruled) for w in check):
        names = {w["pdb"] for w in check}
        rest = [json.loads(json.dumps(w)) for w in work if w["pdb"] not in names]  # leads the checked ones may follow
        apply_rules(check + rest, kb, ruled, *rule_inputs(check + rest, ruled), conf, _conflicts(db, kb, ruled))
    for w in check:
        b = (mem.get(w["pdb"]) or {}).get("barrel")
        if b and "barrelStrands" not in w:
            w["barrelStrands"] = b
    # Draft listings for the whole of each paper: the lead lists every related line, and a
    # related line names what differs from the others.
    drafts = json.loads(json.dumps(work))
    add_listings(drafts, db)
    for w, d in zip(work, drafts):
        new_mp = d.get("mp")
        if not new_mp:
            continue
        merged = {**(w.get("mp") or {}), **new_mp}
        if "relatedLine" not in new_mp:
            merged.pop("relatedLine", None)
        w["mp"] = merged
    by_pdb = {w["pdb"]: w for w in work}
    after = defaultdict(list)
    for w in work:
        after[key_of(w)].append(show(w))
    ordered = cands == sorted(cands, key=lambda c: (-(c.get("membrane") or 0), c["pdb"]))
    doc["candidates"] = [by_pdb.get(c["pdb"], c) for c in cands] + [by_pdb[p] for p in sorted(sibs) if p in by_pdb]
    if ordered:  # the weekly run's order
        doc["candidates"].sort(key=lambda c: (-(c.get("membrane") or 0), c["pdb"]))
    added = [p for p in sorted(sibs) if p in by_pdb]
    stats = doc.setdefault("stats", {})
    stats[siblings_mod.STAT] = stats.get(siblings_mod.STAT, 0) + len(added)
    stats["candidates"] = len(doc["candidates"])
    doc["siblingsAdded"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes")
    path.write_text(json.dumps(doc, separators=(",", ":")), encoding="utf-8")
    out["added"] = [{"pdb": p, "released": by_pdb[p].get("released"), "membrane": by_pdb[p].get("membrane"),
                     "kept": sibs[p]["kept"], "paper": sibs[p]["paper"]} for p in added]
    out["studies"] = [{"paper": k, "before": sorted(before.get(k, [])), "after": sorted(after.get(k, []))}
                      for k in sorted(touched_papers)]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--since", help="first release date to consider (default: 26 weeks ago)")
    r.add_argument("--until")
    sub.add_parser("backtest")
    sub.add_parser("apply-rules", aliases=["replace-rules"])
    sb = sub.add_parser("siblings", help="add same-paper structures released outside the window to data/candidates.json")
    sb.add_argument("--file", help="candidates file to rewrite (default: data/candidates.json)")
    rl = sub.add_parser("relist", help="re-draft the listings in data/candidates.json; nothing is re-scored or re-placed")
    rl.add_argument("--file", help="candidates file to rewrite (default: data/candidates.json)")
    args = ap.parse_args()
    if args.cmd in ("apply-rules", "replace-rules"):
        res = reapply_rules()
        print(json.dumps(res, indent=1))
        return
    db = load(source_xml())
    SITE_DATA.mkdir(exist_ok=True)

    if args.cmd == "relist":
        r = relist(Path(args.file) if args.file else None, db)
        print(f"relist: re-drafted {r['redrafted']} of {r['candidates']} candidates; changed: {r['names']} names, "
              f"{r['descriptions']} descriptions, {r['relatedLines']} related lines")
        return

    if args.cmd == "siblings":
        r = add_siblings(Path(args.file) if args.file else None, db)
        c = r["counts"]
        print(f"siblings: {c.get('papers', 0)} papers of candidates checked in {c.get('seconds', 0)} s; "
              f"{c.get('found', 0)} structures from outside the window not in mpstruc, "
              f"{c.get('belowFloor', 0)} of them under {round(100 * MEMBRANE_FLOOR)}% membrane protein (left out), "
              f"{len(r['added'])} added ({c.get('keptBySameProtein', 0)} kept for having the protein of a candidate "
              f"of the same paper) to {len(r['studies'])} papers")
        for st in r["studies"]:
            print(f"  {st['paper']}\n    before: {'; '.join(st['before'])}\n    after:  {'; '.join(st['after'])}")
        return

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
    # Structures waiting for their paper (data/awaiting.json) are re-checked
    # every week, even after they leave the release window (see awaiting.py).
    from . import autoreject, awaiting
    # Structures curators already approved or rejected (data/curation.json, written
    # by the site's Curators tab) are not suggested again, except rejections
    # waiting for a paper ("Not published yet", see autoreject.py): those stay
    # tracked and checked so they come back to "To decide" when it is out.
    curation_path = SITE_DATA / "curation.json"
    curation = json.loads(curation_path.read_text(encoding="utf-8")) if curation_path.exists() else {}
    decided = autoreject.decided_ids(curation)
    unpub_rejected = set(autoreject.unpublished_rejections(curation))
    wait_path = SITE_DATA / "awaiting.json"
    wait_rec, _ = autoreject.track_unpublished(awaiting.load(wait_path), curation, until)
    wait_rec, given_up = awaiting.drop_stale(wait_rec, until)
    if given_up:
        print(f"awaiting: stopped tracking {len(given_up)} structures still without a paper "
              f"{awaiting.GIVE_UP_DAYS // 365} years after release")
    ids = awaiting.with_tracked(ids, wait_rec)
    try:
        fresh_cits = awaiting.fresh_citations(sorted(awaiting.tracked(wait_rec)))
        awaiting.forget_changed(fresh_cits)
    except Exception as err:  # keep last week's state for these if the PDB is unreachable
        print(f"awaiting: could not re-check tracked structures ({err})")
        fresh_cits = {}
    # Score every new structure for how likely it is to be a membrane protein
    # (mpstruc/membrane.py), and place those at 10% or more.
    from . import membrane
    # mpstruc curators sometimes list a structure only in another entry's
    # description ('alternative conformation, 2.90 Å: 32FD', fragment-screen
    # series), not as an entry or related code; those are already in mpstruc.
    mentioned = mentioned_in_descriptions(db)
    in_text = [p for p in ids if p in mentioned and p not in kb.by_pdb]
    # Every new release that mpstruc already has, with where it is, for the
    # site's 'Already in mpstruc' list and labels.
    already = in_mpstruc(kb, mentioned, ids)
    ids = [p for p in ids if p not in mentioned or p in kb.by_pdb]
    fresh = [p for p in ids if p not in kb.by_pdb and p.upper() not in decided]
    mem = membrane.score_structures(fresh, db)
    floor = MEMBRANE_FLOOR
    keep = {p for p, v in mem.items() if v["p"] >= floor}
    n_scored, n_below = len(mem), sum(v["p"] < floor for v in mem.values())
    # A paper's structures released before or after the window, not in mpstruc and not decided
    # (siblings.py), join the run here, so they are placed with the paper's candidates
    # and tracked, checked against the rules and drafted like them. They are neither in
    # the window's 'released' count nor in the 'Already in mpstruc' list (they are not there).
    sibs, sib_counts = {}, {}
    try:
        # What last week's file knew of each paper: CI keeps no lookup cache between runs,
        # so this is what a paper keeps if the PDB search fails this week.
        try:
            prev = json.loads((SITE_DATA / "candidates.json").read_text(encoding="utf-8")).get("candidates") or []
        except (OSError, ValueError):
            prev = []
        sibs, sib_mem, sib_counts = find_siblings(window_candidates(sorted(keep), fresh_cits), [since, until],
                                                  {p.upper() for p in ids}, kb, mentioned, decided, db, floor,
                                                  fallback=siblings_mod.known_entries(prev))
        mem.update(sib_mem)
        ids = ids + sorted(sibs)
        keep |= set(sibs)
    except Exception as err:  # the weekly run goes on without them
        print(f"siblings: skipped this run ({err})")
        sibs = {}
    props, stats = propose(kb, ids, conf, model, keep=keep)
    stats["released"] -= len(sibs)
    stats["candidates"] = len(props)
    stats["scored"] = n_scored
    stats["already listed"] += len(in_text)
    stats["mentioned only in an mpstruc description"] = len(in_text)
    stats["below 10% membrane protein (not shown)"] = n_below
    stats[siblings_mod.STAT] = sum(p["pdb"] in sibs for p in props)
    props = [_with_membrane(p, mem[p["pdb"]]) for p in props]  # with the per-protein 'Check it yourself' table
    mark_siblings(props, sibs, [since, until])
    props.sort(key=lambda p: (-p["membrane"], p["pdb"]))
    # The curators' rules (curator_rules.json): barrel strands counted in the structure
    # against the typical number of the suggested subfamily.
    ruled = strand_rules(kb)
    conflicts = None
    if ruled:
        conflicts = _conflicts(db, kb, ruled)
        moved = apply_rules(props, kb, ruled, *rule_inputs(props, ruled), conf, conflicts)
        stats["moved by a curators' rule"] = len(moved)
    for p in props:  # a barrel counted for the membrane score is shown on the card too
        if mem[p["pdb"]].get("barrel") and "barrelStrands" not in p:
            p["barrelStrands"] = mem[p["pdb"]]["barrel"]
    # Citations first: any cached entry whose paper turns out to be out now is
    # dropped from the cache, so the draft listing below cites that paper.
    cits = awaiting.citations_for(props, fresh_cits)
    add_listings(props)
    run_day = dt.datetime.now(dt.timezone.utc).date().isoformat()
    wait_rec, wait_stats = awaiting.update(props, {**fresh_cits, **cits}, wait_rec, run_day, hold=unpub_rejected)
    wait_stats["given_up"] = len(given_up)
    awaiting.save(wait_path, wait_rec)
    # Record suggestions still waiting for their paper as rejected by the
    # automation, and bring back those whose paper is now out. Read-modify-
    # write: the site may have saved curator decisions since we read the file.
    cur_text = curation_path.read_text(encoding="utf-8") if curation_path.exists() else None
    cur_new, cur_raw, auto = autoreject.update_text(cur_text, props, wait_rec, autoreject.now_iso())
    if cur_new != cur_text:
        curation_path.write_text(cur_new, encoding="utf-8")
    autoreject.mark_candidates(props, cur_raw)
    autoreject.set_stats(stats, auto)
    stats["awaiting publication"] = wait_stats["awaiting"]
    snap = (db.meta.get("lastDatabaseEditDate") or db.meta.get("lastNewProteinDate") or "")[:10]
    age = (today - dt.date.fromisoformat(snap)).days if snap else None
    out = {"generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
           "window": [since, until], "snapshot": db.meta.get("timeStamp"), "snapshotAgeDays": age, "stats": stats,
           "threshold": 0.5, "confidence": conf, "runDate": run_day,
           "awaiting": wait_stats, "candidates": props, "inMpstruc": already}
    mem_report = ROOT / "reports" / "membrane.json"
    if mem_report.exists():
        mr = json.loads(mem_report.read_text())
        t = mr["test"]
        out["membraneAccuracy"] = {"window": t["window"], "trainWindow": mr["train"]["window"], "n": t["n"],
                                   "auc": t["auc"], "at_50": t["at_50"], "calibration": t["calibration"],
                                   "mpstrucLater": mr.get("mpstruc_later")}
        review = ROOT / "reports" / "membrane_review.json"
        if review.exists():
            rv = json.loads(review.read_text())
            out["membraneAccuracy"]["review"] = {k: rv[k] for k in ("high_score_not_in_databases",
                                                                     "low_score_in_databases", "middle_scores")}
    (SITE_DATA / "candidates.json").write_text(json.dumps(out, separators=(",", ":")))
    try:  # the rules and the official entries that contradict them, for the site and the curators
        write_rule_files(db, kb, out["generated"], conflicts)
    except Exception as err:
        print(f"rules: could not write the rule files ({err})")
    # One line per run for the site's Activity log (kept to the last 200 runs).
    log_path = SITE_DATA / "import_log.json"
    runs = json.loads(log_path.read_text()).get("runs", []) if log_path.exists() else []
    runs.append({"at": out["generated"], "window": [since, until], "released": stats.get("released", 0),
                 "candidates": len(props), "membrane50": sum(p["membrane"] >= 0.5 for p in props),
                 "already_decided": len(decided & set(ids)),
                 "awaiting": wait_stats["awaiting"], "newly_published": wait_stats["newly_published"],
                 "autoRejected": len(auto["rejected"]), "autoUnrejected": len(auto["unrejected"]),
                 "siblings": stats[siblings_mod.STAT], "siblingSeconds": sib_counts.get("seconds", 0)})
    log_path.write_text(json.dumps({"runs": runs[-200:]}, indent=1))
    print(f"{stats}; {sum(p['membrane'] >= 0.5 for p in props)} at least 50% likely to be membrane proteins; "
          f"wrote {SITE_DATA / 'candidates.json'}")


if __name__ == "__main__":
    main()
