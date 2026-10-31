"""How likely is a PDB structure to be a membrane protein?

    python -m mpstruc.membrane        # train, test, write membrane_model.json + reports/membrane.md

Labels (what counts as a membrane protein): the structure is listed in
mpstruc, or PDBTM or MemProtMD (which screened every PDB entry for
transmembrane proteins) list it, or OPM classes it as transmembrane.
Structures OPM classes only as monotopic/peripheral or peptide, and that no
other source lists, are left out of training and testing (OPM does not say
whether they sit in the membrane or on it). Everything else released in the
same weeks is labelled "not a membrane protein". PDBTM screened the whole PDB
consistently until March 2025, so training (Jul 2023 - Aug 2024) and testing
(Sep 2024 - Mar 2025) use that period.

Features are things known on the day a structure is released, so the
prediction works for brand-new entries: whether the residues that cross the
membrane (UniProt transmembrane segments) are actually resolved in the
structure, UniProt keywords and families,
detergents, lipids or nanodiscs in the sample, whether mpstruc already lists
the same protein (only earlier mpstruc entries count), how often its protein
family is a membrane protein in the training data, the method, and words in
the PDB title. The probability is calibrated: of structures given about 80%,
about 80% are membrane proteins in the test period (see reports/membrane.md).
"""
from __future__ import annotations

import json
import math
import re
import urllib.request
from collections import Counter, defaultdict

from . import membranedb, rcsb, uniprot
from .automation import FUSION_UNIPROT, SEARCH, TM_KEYWORDS, released_between
from .model import load
from .report import ROOT, source_xml

TRAIN = ("2023-07-01", "2024-08-31")
TEST = ("2024-09-01", "2025-03-31")
MODEL_PATH = ROOT / "membrane_model.json"

DETERGENTS = set("""LMT LMN BOG C8E LDA UMQ DMU BNG HTG CPS SDS UDM DDQ J2K AV0 PX4 OLC 78M 78N 7E9 MPG
Y01 Q7G 2J8 UNL LDM GDN D12 TRD""".split()) - {"UNL"}
LIPIDS = set("""PCW POV PEE PLC LPP PGW CDL 3PE PC1 6PL P5S PEF PEV PTY LHG LMG DGD SQD CLR PA8 8PE 44E
2CV PCF PGT 7PH PTD 17F J40 P0E LBN 9Z9 PGV PSF""".split())
NANODISC = re.compile(r"membrane scaffold|apolipoprotein a-?i\b|saposin|msp1|msp2", re.I)
# UniProt keywords that point away from (or towards) the membrane.
KW_SOLUBLE = {"KW-0964": "secreted", "KW-0963": "in the cytoplasm", "KW-0539": "in the nucleus"}
KW_MEMBRANE = {"KW-0472", "KW-1003", "KW-0999", "KW-1000", "KW-0256", "KW-0333", "KW-0998", "KW-0793"}
TITLE_STOP = {"of", "the", "in", "with", "and", "a", "an", "from", "to", "by", "for", "on", "at", "bound",
              "structure", "structures", "complex", "state", "form", "cryo-em", "cryoem", "crystal", "x-ray"}


def _search(q: dict) -> set[str]:
    req = urllib.request.Request(SEARCH, data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        body = r.read()
    return {x["identifier"] for x in json.loads(body).get("result_set", [])} if body else set()


def annotated(start: str, end: str, types: list[str]) -> set[str]:
    """PDB entries released in [start, end] that the PDB tags with any of these annotation types."""
    return _search({"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_accession_info.initial_release_date", "operator": "range",
            "value": {"from": start, "to": end, "include_lower": True, "include_upper": True}}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_annotation.type", "operator": "in", "value": types}}]},
        "return_type": "entry", "request_options": {"return_all_hits": True}})


# ---------------------------------------------------------------- features

class Context:
    """What mpstruc and the training data know, for computing features."""

    def __init__(self, db, fam_stats: dict[str, list[int]] | None = None):
        self.mpstruc_ids = {e.pdbCode.upper() for e in db.entries()} | {r.pdb.upper() for e in db.entries() for r in e.related}
        ids = sorted(self.mpstruc_ids)
        ext, det = rcsb.entry_extras(ids), rcsb.entry_entities(ids)  # entities carry the UniProt IDs
        # UniProt accession -> release dates of mpstruc structures containing it.
        self.protein_dates: dict[str, list[str]] = defaultdict(list)
        self.protein_example: dict[str, str] = {}
        for pdb in ids:
            x, d = ext.get(pdb) or {}, det.get(pdb)
            rel = x.get("released") or ""
            for u in _uniprots(d) - FUSION_UNIPROT:
                self.protein_dates[u].append(rel)
                self.protein_example.setdefault(u, pdb)
        for v in self.protein_dates.values():
            v.sort()
        self.fam_stats = fam_stats or {}  # InterPro -> [membrane, not membrane] counts in training
        self.tcdb = membranedb.tcdb()


def _uniprots(det: dict | None) -> set[str]:
    return {u for e in (det or {}).get("polymer_entities") or []
            for u in (e.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []}


def _fetch(ids: list[str]):
    ext, fam, det, cov = rcsb.entry_extras(ids), rcsb.entry_families(ids), rcsb.entry_details(ids), rcsb.entry_coverage(ids)
    accs = sorted({a for v in cov.values() for e in (v or []) for a in e.get("align") or {}})
    return ext, fam, det, cov, uniprot.tm_segments(accs)


def tm_coverage(cov: list | None, tmseg: dict[str, list]) -> tuple[int, int]:
    """(membrane-crossing segments in the construct, of those actually modelled).

    Maps each UniProt transmembrane/intramembrane segment onto the chain via the
    PDB's UniProt alignment, and counts it as modelled if at least half of its
    residues are not in the chain's unobserved ranges.
    """
    total = modelled = 0
    for ent in cov or []:
        unobs = ent.get("unobs") or []
        for acc, regions in (ent.get("align") or {}).items():
            for b, e in tmseg.get(acc) or []:
                pos = []
                for u in range(b, e + 1):
                    for eb, rb, n in regions:
                        if rb <= u < rb + n:
                            pos.append(eb + u - rb)
                            break
                if len(pos) < (e - b + 1) / 2:
                    continue  # this part of the protein is not in the construct
                total += 1
                seen = sum(1 for q in pos if not any(lo <= q <= hi for lo, hi in unobs))
                modelled += seen >= (e - b + 1) / 2
    return total, modelled


def raw_facts(pdb: str, ext: dict | None, fam: dict | None, det: dict | None,
              cov: list | None = None, tmseg: dict | None = None) -> dict | None:
    """Release-day facts about a structure (no labels, no PDB annotations)."""
    if not det:
        return None
    kw = (ext or {}).get("kw") or {}
    prot = [e for e in det.get("polymer_entities") or []
            if (e.get("entity_poly") or {}).get("rcsb_entity_polymer_type") == "Protein"]
    ups = _uniprots(det) - FUSION_UNIPROT
    tm_ups = {u for u in ups if TM_KEYWORDS & set(kw.get(u, ()))}
    all_kw = {k for u in ups for k in kw.get(u, ())}
    ligs = {((n.get("nonpolymer_comp") or {}).get("chem_comp") or {}).get("id") for n in det.get("nonpolymer_entities") or []}
    descs = [((e.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "") for e in prot]
    ipr = (fam or {}).get("ipr") or {}
    has_tm_protein = any(tmseg.get(u) for u in ups) if tmseg is not None else None
    tm_total, tm_modelled = tm_coverage(cov, tmseg or {}) if cov is not None else (None, None)
    return {
        "has_tm_protein": has_tm_protein, "tm_in_construct": tm_total, "tm_modelled": tm_modelled,
        "pdb": pdb, "released": (ext or {}).get("released", ""),
        "uniprots": sorted(ups), "tm_uniprots": sorted(tm_ups), "keywords": sorted(all_kw),
        "detergents": sorted(ligs & DETERGENTS), "lipids": sorted(ligs & LIPIDS),
        "nanodisc": any(NANODISC.search(d) for d in descs),
        "ipr": sorted({x for u in ups for x in ipr.get(u, ())}),
        "methods": " ".join(m.get("method", "") for m in det.get("exptl") or []).upper(),
        "title": ((det.get("struct") or {}).get("title") or ""),
        "n_protein": len(prot),
        "n_unmapped": sum(1 for e in prot if not (e.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids")),
    }


def _fam_logodds(stats: list[int] | None, prior: float) -> float:
    """Smoothed log-odds that a structure with this InterPro family is a membrane protein."""
    pos, neg = stats or (0, 0)
    a = 2.0  # pseudo-counts at the base rate
    p = (pos + a * prior) / (pos + neg + a)
    return math.log(p / (1 - p)) - math.log(prior / (1 - prior))


def features(f: dict, ctx: Context, prior: float, own_label: int | None = None) -> tuple[dict[str, float], dict]:
    """Model features plus the facts used to explain them.

    own_label: for training rows, the structure's own label is removed from
    the family counts, so a structure never votes for itself.
    """
    x: dict[str, float] = {}
    why: dict = {}
    kw = set(f["keywords"])
    x["kw_tm"] = float(bool(f["tm_uniprots"]))
    x["kw_tm_all"] = float(bool(f["uniprots"]) and len(f["tm_uniprots"]) == len(f["uniprots"]))
    x["kw_membrane"] = float(bool(kw & KW_MEMBRANE))
    for k in KW_SOLUBLE:
        x[f"kw_{k}"] = float(k in kw)
    x["no_uniprot"] = float(not f["uniprots"])
    # Is the membrane-crossing part actually in the structure? (None: coverage unknown.)
    if f.get("tm_modelled") is not None:
        x["tm_modelled"] = float(f["tm_modelled"] > 0)
        x["tm_missing"] = float(bool(f.get("has_tm_protein")) and f["tm_modelled"] == 0)
    tcs = sorted({tc for u in f["uniprots"] for tc in ctx.tcdb.get(u, ())})
    x["tcdb_transporter"] = float(any(tc[0] in "123" for tc in tcs))
    x["tcdb_other"] = float(bool(tcs) and not x["tcdb_transporter"])
    why["tcdb"] = tcs[:2]
    x["detergent"] = float(bool(f["detergents"]))
    x["lipid"] = float(bool(f["lipids"]))
    x["nanodisc"] = float(f["nanodisc"])
    earlier = [u for u in f["uniprots"] if any(d and d < f["released"] for d in ctx.protein_dates.get(u, ()))]
    x["same_protein"] = float(bool(earlier))
    why["same_protein"] = ctx.protein_example.get(earlier[0]) if earlier else None
    best, best_ipr = 0.0, None
    for i in f["ipr"]:
        st = list(ctx.fam_stats.get(i, (0, 0)))
        if own_label is not None and sum(st):
            st[0 if own_label else 1] = max(0, st[0 if own_label else 1] - 1)
        lo = _fam_logodds(st, prior)
        if abs(lo) > abs(best):
            best, best_ipr = lo, (i, st)
    x["family"] = max(-4.0, min(4.0, best))
    why["family"] = best_ipr
    m = f["methods"]
    x["m:em"], x["m:xray"], x["m:nmr"] = float("ELECTRON MICROSCOPY" in m), float("X-RAY" in m), float("NMR" in m)
    x["many_chains"] = float(f["n_protein"] > 4)
    for w in set(re.findall(r"[a-z][a-z0-9-]{2,}", f["title"].lower())) - TITLE_STOP:
        x[f"t:{w}"] = 1.0
    return x, why


# ---------------------------------------------------------------- model

def fit(X: list[dict], y: list[int], l2: float = 10.0, steps: int = 3000, min_count: int = 40) -> dict:
    """Plain (unweighted) logistic regression, so the output is a probability."""
    import numpy as np
    counts = Counter(k for x in X for k, v in x.items() if v)
    keys = sorted(k for k in counts if not k.startswith("t:") or counts[k] >= min_count)
    A = np.array([[x.get(k, 0.0) for k in keys] for x in X], dtype=float)
    t = np.array(y, dtype=float)
    w, b = np.zeros(len(keys)), math.log(t.mean() / (1 - t.mean()))
    for _ in range(steps):
        p = 1 / (1 + np.exp(-(A @ w + b)))
        g = p - t
        w -= 0.5 * ((A.T @ g) / len(t) + l2 * w / len(t))
        b -= 0.5 * g.mean()
    return {"bias": float(b), "weights": {k: float(v) for k, v in zip(keys, w) if abs(v) >= 0.005}}


def predict(x: dict[str, float], model: dict) -> float:
    z = model["bias"] + sum(model["weights"].get(k, 0.0) * v for k, v in x.items())
    p = 1 / (1 + math.exp(-z))
    if model.get("isotonic"):  # observed frequencies on held-out training structures
        return _apply_isotonic(model["isotonic"], p)
    cal = model.get("calibration")
    if cal:  # Platt scaling fitted on held-out weeks
        z2 = cal["a"] * math.log(max(p, 1e-6) / max(1 - p, 1e-6)) + cal["b"]
        p = 1 / (1 + math.exp(-z2))
    return p


def explain(x: dict[str, float], why: dict, f: dict, p: float, model: dict, n: int = 3) -> list[str]:
    """Plain-English reasons, strongest first: for a high score the evidence
    for, then the strongest against; for a low score the other way round."""
    w = model["weights"]
    reasons = []  # (signed contribution, text)
    def add(key, val, text):
        if text:
            reasons.append((w.get(key, 0.0) * val, text))
    kw = set(f["keywords"])
    if x.get("kw_tm"):
        add("kw_tm", 1, "UniProt lists " + ("all its proteins" if x.get("kw_tm_all") and len(f["uniprots"]) > 1
                                            else "a protein in it") + " as crossing the membrane")
    elif f["uniprots"]:
        reasons.append((-abs(w.get("kw_tm", 1.0)), "UniProt does not list any of its proteins as crossing the membrane"))
    if x.get("kw_membrane") and not x.get("kw_tm"):
        add("kw_membrane", 1, "UniProt places a protein in it in a membrane")
    for k, place in KW_SOLUBLE.items():
        if x.get(f"kw_{k}"):
            add(f"kw_{k}", 1, f"UniProt says a protein in it is {place}")
    if x.get("tcdb_transporter"):
        add("tcdb_transporter", 1, f"TCDB classifies a protein in it as a channel or transporter (TC {why['tcdb'][0]})")
    elif x.get("tcdb_other"):
        add("tcdb_other", 1, f"TCDB lists a protein in it (TC {why['tcdb'][0]})")
    if x.get("tm_modelled"):
        add("tm_modelled", 1, f"{f['tm_modelled']} of its membrane-crossing segments are resolved in the structure")
    if x.get("tm_missing"):
        add("tm_missing", 1, "its membrane-crossing segments are not resolved, so only an outside part of the protein is in the structure")
    if x.get("detergent"):
        add("detergent", 1, f"the sample contains detergent ({', '.join(f['detergents'][:2])})")
    if x.get("lipid"):
        add("lipid", 1, f"lipids are bound ({', '.join(f['lipids'][:2])})")
    if x.get("nanodisc"):
        add("nanodisc", 1, "it was solved in a nanodisc")
    if x.get("same_protein"):
        add("same_protein", 1, f"mpstruc already lists this protein (e.g. {why['same_protein']})")
    fam = why.get("family")
    if fam and abs(x.get("family", 0)) > 0.5:
        i, (pos, neg) = fam
        add("family", x["family"], f"its protein family ({i}) is a membrane protein in {pos} of {pos + neg} training structures")
    words = [(w.get(k, 0.0), k[2:]) for k in x if k.startswith("t:") and abs(w.get(k, 0.0)) > 0.3]
    for wt, word in sorted(words, key=lambda t: -abs(t[0]))[:1]:
        reasons.append((wt, f"the title mentions '{word}'"))
    if not f["uniprots"]:
        reasons.append((w.get("no_uniprot", 0.0), "none of its proteins is linked to UniProt yet"))
    high = p >= 0.5
    pro = sorted([r for r in reasons if (r[0] > 0) == high and abs(r[0]) > 0.05], key=lambda r: -abs(r[0]))
    con = sorted([r for r in reasons if (r[0] > 0) != high and abs(r[0]) > 0.05], key=lambda r: -abs(r[0]))
    out = [t for _, t in pro[:n]]
    if con:
        out.append("but " + con[0][1])
    return out


def _labels(start: str, end: str, mpstruc_ids: set[str]) -> tuple[list[str], set[str]]:
    """Release IDs (minus OPM-only monotopic/peripheral and peptide entries) and the membrane-protein set."""
    ids = released_between(start, end)
    opm = membranedb.opm()
    pos = annotated(start, end, ["PDBTM", "MemProtMD"]) | (set(ids) & mpstruc_ids) | \
        {p for p in ids if opm.get(p) == "transmembrane"}
    unsure = {p for p in ids if opm.get(p) in ("monotopic/peripheral", "peptide")} - pos
    return [p for p in ids if p not in unsure], pos


def evaluate(rows, model) -> dict:
    return evaluate_preds([(predict(x, model), y) for _, x, y, _, _ in rows])


def evaluate_preds(preds: list[tuple[float, int]]) -> dict:
    n = len(preds)
    bins = []
    for lo in [i / 10 for i in range(10)]:
        b = [(p, y) for p, y in preds if lo <= p < lo + 0.1 or (lo == 0.9 and p == 1.0)]
        if b:
            bins.append({"range": f"{int(lo * 100)}-{int(lo * 100) + 10}%", "n": len(b),
                         "predicted": round(100 * sum(p for p, _ in b) / len(b), 1),
                         "actual": round(100 * sum(y for _, y in b) / len(b), 1)})
    pos = sum(y for _, y in preds)
    tp = sum(1 for p, y in preds if p >= 0.5 and y)
    fp = sum(1 for p, y in preds if p >= 0.5 and not y)
    ranked = sorted(preds, key=lambda t: -t[0])
    # AUC by ranking (Mann-Whitney)
    rank_sum, r = 0.0, 0
    for i, (p, y) in enumerate(sorted(preds, key=lambda t: t[0]), 1):
        if y:
            rank_sum += i
    auc = (rank_sum - pos * (pos + 1) / 2) / (pos * (n - pos)) if pos and n - pos else None
    brier = sum((p - y) ** 2 for p, y in preds) / n
    return {"n": n, "membrane": pos, "auc": round(auc, 4) if auc else None, "brier": round(brier, 4),
            "at_50": {"called_membrane": tp + fp, "precision_pct": round(100 * tp / (tp + fp), 1) if tp + fp else None,
                      "recall_pct": round(100 * tp / pos, 1) if pos else None},
            "calibration": bins, "top_50_precision": round(100 * sum(y for _, y in ranked[:50]) / 50, 1)}


def _platt(preds: list[tuple[float, int]]) -> dict:
    """Fit p' = sigmoid(a * logit(p) + b) on held-out predictions."""
    import numpy as np
    z = np.array([math.log(max(p, 1e-6) / max(1 - p, 1e-6)) for p, _ in preds])
    t = np.array([y for _, y in preds], dtype=float)
    a, b = 1.0, 0.0
    for _ in range(5000):
        q = 1 / (1 + np.exp(-(a * z + b)))
        a -= 0.05 * float(((q - t) * z).mean())
        b -= 0.05 * float((q - t).mean())
    return {"a": a, "b": b}


def _isotonic(preds: list[tuple[float, int]], bins: int = 40) -> list[list[float]]:
    """Pool-adjacent-violators fit of P(membrane | raw score) on held-out
    predictions. Returns [[raw score, calibrated probability], ...] knots."""
    pts = sorted(preds)
    size = max(1, len(pts) // bins)
    prior = sum(y for _, y in preds) / len(preds)
    blocks = []  # [sum_y, n, mean_x]; one pseudo-observation at the base rate keeps rates inside (0, 1)
    for i in range(0, len(pts), size):
        chunk = pts[i:i + size]
        if blocks and len(chunk) < size // 2:  # fold a small remainder into the last block
            y0, n0, x0 = blocks.pop()
            chunk_y, chunk_n = sum(y for _, y in chunk), len(chunk)
            blocks.append([y0 + chunk_y, n0 + chunk_n, (x0 * (n0 - 1) + sum(p for p, _ in chunk)) / (n0 - 1 + chunk_n)])
            continue
        blocks.append([sum(y for _, y in chunk) + prior, len(chunk) + 1, sum(p for p, _ in chunk) / len(chunk)])
    merged: list[list[float]] = []
    for b in blocks:
        merged.append(b[:])
        while len(merged) > 1 and merged[-2][0] / merged[-2][1] > merged[-1][0] / merged[-1][1]:
            y2, n2, x2 = merged.pop()
            y1, n1, x1 = merged.pop()
            merged.append([y1 + y2, n1 + n2, (x1 * n1 + x2 * n2) / (n1 + n2)])
    return [[x, y / n] for y, n, x in merged]


def _apply_isotonic(knots: list[list[float]], p: float) -> float:
    if p <= knots[0][0]:
        return knots[0][1]
    for (x0, y0), (x1, y1) in zip(knots, knots[1:]):
        if p <= x1:
            return y0 + (y1 - y0) * (p - x0) / (x1 - x0) if x1 > x0 else y1
    return knots[-1][1]


def _family_stats(rows) -> dict[str, list[int]]:
    stats: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for f, y in rows:
        for i in set(f["ipr"]):
            stats[i][0 if y else 1] += 1
    return dict(stats)


def train_and_test(folds: int = 5) -> dict:
    """Train on TRAIN, calibrate on held-out folds of it, report on TEST.

    Family statistics are computed without the structure being scored (from
    the other time blocks of the training period), and the calibration maps
    out-of-fold predictions to observed frequencies, so the percentage shown
    is what was actually seen on structures the model had not learnt from.
    """
    db = load(source_xml())
    ctx = Context(db)
    ids_tr, pos_tr = _labels(*TRAIN, ctx.mpstruc_ids)
    ext, fam, det, cov, seg = _fetch(ids_tr)
    facts = []
    for pdb in ids_tr:
        f = raw_facts(pdb, ext.get(pdb), fam.get(pdb), det.get(pdb), cov.get(pdb), seg)
        if f:
            facts.append((f, int(pdb in pos_tr)))
    facts.sort(key=lambda r: (r[0]["released"], r[0]["pdb"]))
    prior = sum(y for _, y in facts) / len(facts)
    block = [min(folds - 1, i * folds // len(facts)) for i in range(len(facts))]
    # Features for each training structure use family statistics from the other blocks only.
    rows = []
    for k in range(folds):
        ctx.fam_stats = _family_stats([r for r, b in zip(facts, block) if b != k])
        rows += [(f["pdb"], *features(f, ctx, prior), y, f, k) for (f, y), b in zip(facts, block) if b == k]
    # rows: (pdb, x, why, y, facts, block)
    oof = []
    for k in range(folds):
        m = fit([r[1] for r in rows if r[5] != k], [r[3] for r in rows if r[5] != k])
        oof += [(predict(r[1], m), r[3]) for r in rows if r[5] == k]
    model = fit([r[1] for r in rows], [r[3] for r in rows])
    model["isotonic"] = _isotonic(oof)
    model["prior"] = prior
    model["fam_stats"] = {k: v for k, v in _family_stats(facts).items() if sum(v) >= 2}
    ctx.fam_stats = model["fam_stats"]

    ids_te, pos_te = _labels(*TEST, ctx.mpstruc_ids)
    ext, fam, det, cov, seg = _fetch(ids_te)
    rows_te = []
    for pdb in ids_te:
        f = raw_facts(pdb, ext.get(pdb), fam.get(pdb), det.get(pdb), cov.get(pdb), seg)
        if f:
            x, why = features(f, ctx, prior)
            rows_te.append((pdb, x, int(pdb in pos_te), f, why))
    oof_eval = [(_apply_isotonic(model["isotonic"], p), y) for p, y in oof]
    result = {"train": {"window": TRAIN, **evaluate_preds(oof_eval)},
              "test": {"window": TEST, **evaluate(rows_te, model)}}
    MODEL_PATH.write_text(json.dumps(model))
    later = sorted(set(released_between("2025-07-01", "2025-12-31")) & ctx.mpstruc_ids)
    sc = score_structures(later, db, model)
    ps = [v["p"] for v in sc.values()]
    result["mpstruc_later"] = {"window": ["2025-07-01", "2025-12-31"], "n": len(ps),
                               "at_least_50_pct": round(100 * sum(p >= 0.5 for p in ps) / len(ps), 1) if ps else None,
                               "at_least_10_pct": round(100 * sum(p >= 0.1 for p in ps) / len(ps), 1) if ps else None}
    examples = sorted(rows_te, key=lambda r: -predict(r[1], model))
    result["examples"] = [{"pdb": r[0], "membrane": r[2], "p": round(predict(r[1], model), 3),
                           "why": explain(r[1], r[4], r[3], predict(r[1], model), model)} for r in examples[::max(1, len(examples) // 40)]]
    (ROOT / "reports" / "membrane.json").write_text(json.dumps(result, indent=1))
    (ROOT / "reports" / "membrane.md").write_text(report_md(result, model))
    return result


def load_model() -> dict | None:
    return json.loads(MODEL_PATH.read_text()) if MODEL_PATH.exists() else None


def score_structures(ids: list[str], db=None, model: dict | None = None) -> dict[str, dict]:
    """{pdb: {"p": probability, "why": [reasons]}} for new structures."""
    model = model or load_model()
    if not model:
        return {}
    ctx = Context(db or load(source_xml()), model.get("fam_stats"))
    ext, fam, det, cov, seg = _fetch(ids)
    out = {}
    for pdb in ids:
        f = raw_facts(pdb, ext.get(pdb), fam.get(pdb), det.get(pdb), cov.get(pdb), seg)
        if not f:
            continue
        x, why = features(f, ctx, model["prior"])
        p = predict(x, model)
        out[pdb] = {"p": round(p, 3), "why": explain(x, why, f, p, model)}
    return out


def report_md(r: dict, model: dict) -> str:
    t = r["test"]
    L = ["# How likely is a structure to be a membrane protein?", "",
         "A structure counts as a membrane protein if mpstruc lists it, if PDBTM or MemProtMD (which screened every "
         "PDB entry for transmembrane proteins) do, or if OPM classes it as transmembrane. Structures OPM classes only as "
         "monotopic/peripheral or peptide, with no other source, are left out. Signals include UniProt keywords and "
         "families and TCDB transporter classes. The model uses only facts known on the day a structure is released. "
         f"It was trained on PDB releases {TRAIN[0]} to {TRAIN[1]} (calibrated on structures each fold had not seen) and "
         f"tested on {TEST[0]} to {TEST[1]}.", "",
         "## Test period", "",
         f"- Structures: {t['n']:,}, of which membrane proteins: {t['membrane']:,}",
         f"- Ranking quality (AUC): {t['auc']}; Brier score: {t['brier']}",
         f"- Called membrane protein (≥50%): {t['at_50']['called_membrane']:,}; of those, {t['at_50']['precision_pct']}% are; "
         f"they include {t['at_50']['recall_pct']}% of all membrane proteins",
         f"- Of the {r['mpstruc_later']['n']} structures mpstruc added from Jul–Dec 2025 releases, "
         f"{r['mpstruc_later']['at_least_50_pct']}% score at least 50% and {r['mpstruc_later']['at_least_10_pct']}% at least 10%",
         "", "### Does the percentage mean what it says?", "",
         "| Shown as | Structures | Average shown | Actually membrane proteins |", "|---|---|---|---|"]
    L += [f"| {b['range']} | {b['n']:,} | {b['predicted']}% | {b['actual']}% |" for b in t["calibration"]]
    L += ["", "## Strongest signals", "", "| Signal | Weight |", "|---|---|"]
    named = sorted(((k, v) for k, v in model["weights"].items() if not k.startswith("t:")), key=lambda kv: -abs(kv[1]))
    L += [f"| {k} | {v:+.2f} |" for k, v in named]
    words = sorted(((k[2:], v) for k, v in model["weights"].items() if k.startswith("t:")), key=lambda kv: -kv[1])
    L += ["", "Title words that raise it: " + ", ".join(f"{k} ({v:+.2f})" for k, v in words[:15]),
          "", "Title words that lower it: " + ", ".join(f"{k} ({v:+.2f})" for k, v in words[::-1][:15]), "",
          "## Examples from the test period", "", "| PDB | Shown | Membrane protein? | Why |", "|---|---|---|---|"]
    L += [f"| {e['pdb']} | {round(100 * e['p'])}% | {'yes' if e['membrane'] else 'no'} | {'; '.join(e['why'])} |" for e in r["examples"]]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    res = train_and_test()
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != 'calibration'} for k, v in res.items() if k != "examples"}, indent=1))
    print(json.dumps(res["test"]["calibration"], indent=1))
