"""How likely is a PDB structure to be a membrane protein?

    python -m mpstruc.membrane        # train, test, write membrane_model.json + reports/membrane.md

Labels (what counts as a membrane protein): the structure is listed in
mpstruc, or the PDBTM or MemProtMD databases, which screened every PDB entry
for transmembrane proteins, list it. Everything else released in the same
weeks is labelled "not a membrane protein". PDBTM screened the whole PDB
consistently until March 2025, so training (Jul 2023 - Aug 2024) and testing
(Sep 2024 - Mar 2025) use that period.

Features are things known on the day a structure is released, so the
prediction works for brand-new entries: UniProt keywords and families,
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

from . import rcsb
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


def _uniprots(det: dict | None) -> set[str]:
    return {u for e in (det or {}).get("polymer_entities") or []
            for u in (e.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []}


def raw_facts(pdb: str, ext: dict | None, fam: dict | None, det: dict | None) -> dict | None:
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
    return {
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

def fit(X: list[dict], y: list[int], l2: float = 3.0, steps: int = 3000, min_count: int = 20) -> dict:
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
    ids = released_between(start, end)
    pos = annotated(start, end, ["PDBTM", "MemProtMD"]) | (set(ids) & mpstruc_ids)
    return ids, pos


def _rows(ids, pos, ctx, prior, training: bool):
    ext, fam, det = rcsb.entry_extras(ids), rcsb.entry_families(ids), rcsb.entry_details(ids)
    rows = []
    for pdb in ids:
        f = raw_facts(pdb, ext.get(pdb), fam.get(pdb), det.get(pdb))
        if not f:
            continue
        y = int(pdb in pos)
        x, why = features(f, ctx, prior, own_label=y if training else None)
        rows.append((pdb, x, y, f, why))
    return rows


def evaluate(rows, model) -> dict:
    preds = [(predict(x, model), y) for _, x, y, _, _ in rows]
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


def train_and_test() -> dict:
    db = load(source_xml())
    ctx = Context(db)
    ids_tr, pos_tr = _labels(*TRAIN, ctx.mpstruc_ids)
    prior = len(pos_tr & set(ids_tr)) / len(ids_tr)
    # Family statistics from the training period only.
    fam = rcsb.entry_families(ids_tr)
    stats: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for pdb in ids_tr:
        iprs = {x for v in ((fam.get(pdb) or {}).get("ipr") or {}).values() for x in v}
        for i in iprs:
            stats[i][0 if pdb in pos_tr else 1] += 1
    ctx.fam_stats = dict(stats)
    rows_tr = _rows(ids_tr, pos_tr, ctx, prior, training=True)
    # Fit on all but the last two training months, calibrate on those, then refit on everything.
    cut = "2024-07-01"
    early = [r for r in rows_tr if r[3]["released"] < cut]
    late = [r for r in rows_tr if r[3]["released"] >= cut]
    m0 = fit([r[1] for r in early], [r[2] for r in early])
    cal = _platt([(predict(r[1], m0), r[2]) for r in late])
    model = fit([r[1] for r in rows_tr], [r[2] for r in rows_tr])
    model["calibration"] = cal
    model["prior"] = prior
    ids_te, pos_te = _labels(*TEST, ctx.mpstruc_ids)
    rows_te = _rows(ids_te, pos_te, ctx, prior, training=False)
    result = {"train": {"window": TRAIN, **evaluate(rows_tr, model)},
              "test": {"window": TEST, **evaluate(rows_te, model)}}
    model["fam_stats"] = {k: v for k, v in stats.items() if sum(v) >= 2}
    MODEL_PATH.write_text(json.dumps(model))
    # How the score treats what mpstruc actually added later (Jul-Dec 2025, after the PDBTM period).
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
    ext, fam, det = rcsb.entry_extras(ids), rcsb.entry_families(ids), rcsb.entry_details(ids)
    out = {}
    for pdb in ids:
        f = raw_facts(pdb, ext.get(pdb), fam.get(pdb), det.get(pdb))
        if not f:
            continue
        x, why = features(f, ctx, model["prior"])
        p = predict(x, model)
        out[pdb] = {"p": round(p, 3), "why": explain(x, why, f, p, model)}
    return out


def report_md(r: dict, model: dict) -> str:
    t = r["test"]
    L = ["# How likely is a structure to be a membrane protein?", "",
         "A structure counts as a membrane protein if mpstruc lists it, or if PDBTM or MemProtMD (which screened every "
         "PDB entry for transmembrane proteins) do. The model uses only facts known on the day a structure is released. "
         f"It was trained on PDB releases {TRAIN[0]} to {TRAIN[1]} and tested on {TEST[0]} to {TEST[1]}.", "",
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
