"""How likely is a PDB structure to be a membrane protein?

    python -m mpstruc.membrane             # train, test, write membrane_model.json + reports/membrane.md
    python -m mpstruc.membrane --tradeoff  # also compare ways of labelling monotopic proteins
    python -m mpstruc.membrane table       # add the per-protein evidence table (tmEvidence) to data/candidates.json

Labels (what counts as a membrane protein): the structure is listed in
mpstruc, or PDBTM or MemProtMD (which screened every PDB entry for
transmembrane proteins) list it, or OPM classes it as transmembrane.
Structures mpstruc lists only as monotopic are handled by MONOTOPIC_LABEL
(see _labels; --tradeoff measures the alternatives). Structures OPM classes
only as monotopic/peripheral or peptide, and that no other source lists, are
left out of training and testing. Everything else released in the same weeks
is labelled "not a membrane protein". PDBTM screened the whole PDB
consistently until March 2025, so training (Jul 2023 - Aug 2024) and testing
(Sep 2024 - Mar 2025) use that period.

Features are things known on the day a structure is released, so the
prediction works for brand-new entries. The strongest is whether the part
that crosses the membrane is resolved in the structure: UniProt transmembrane
segments (helices, or strands where UniProt maps them), or, for beta-barrels
UniProt maps no strands for, a UniProt barrel domain such as 'TBDR
beta-barrel' (uniprot.membrane_features, barrel_domains). A structure of only
a soluble part (a plug, an ectodomain) does not count. Also: whether UniProt
calls a protein transmembrane but maps nothing to check, UniProt keywords
(including "Transmembrane beta strand" and "Cell outer membrane") and
families, TCDB classes (including 1.B, outer-membrane beta-barrel porins),
detergents, lipids or nanodiscs in the sample, whether mpstruc already lists
the same protein (only earlier mpstruc entries count), how often its protein
family is a membrane protein in the training data, the method, and words in
the PDB title. UniProt annotations are read as they are today, as before.
The probability is calibrated: of structures given about 80%, about 80% are
membrane proteins in the test period (see reports/membrane.md).

The curators' rule: a structure is a membrane protein when a significant
membrane-crossing part is resolved in it. When CORE_MIN_SEGMENTS or more UniProt
transmembrane segments (at least CORE_MIN_LENGTH residues) of a main protein are
resolved, or a membrane beta-barrel is, the score is at least TM_FLOOR, so family
statistics or 'secreted' keywords of a partner cannot pull it below 50%, and the
resolved part is the first reason given (core_tm, TM_FLOOR).
"""
from __future__ import annotations

import datetime
import json
import difflib
import functools
import math
import re
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from . import membranedb, rcsb, uniprot
from .automation import FUSION_UNIPROT, SEARCH, TM_KEYWORDS, released_between
from .model import load
from .report import ROOT, source_xml

TRAIN = ("2023-07-01", "2024-08-31")
# How structures mpstruc lists only as monotopic are labelled ("positive", "exclude" or
# "negative"; see _labels and the trade-off table in reports/membrane.md).
MONOTOPIC_LABEL = "positive"
TEST = ("2024-09-01", "2025-03-31")
MODEL_PATH = ROOT / "membrane_model.json"

DETERGENTS = set("""LMT LMN BOG C8E LDA UMQ DMU BNG HTG CPS SDS UDM DDQ J2K AV0 PX4 OLC 78M 78N 7E9 MPG
Y01 Q7G 2J8 UNL LDM GDN D12 TRD""".split()) - {"UNL"}
LIPIDS = set("""PCW POV PEE PLC LPP PGW CDL 3PE PC1 6PL P5S PEF PEV PTY LHG LMG DGD SQD CLR PA8 8PE 44E
2CV PCF PGT 7PH PTD 17F J40 P0E LBN 9Z9 PGV PSF""".split())
NANODISC = re.compile(r"membrane scaffold|apolipoprotein a-?i\b|saposin|msp1|msp2", re.I)
# UniProt keywords that point away from (or towards) the membrane.
KW_SOLUBLE = {"KW-0964": "secreted", "KW-0963": "in the cytoplasm", "KW-0539": "in the nucleus"}
KW_STRAND, KW_OUTER, KW_SECRETED = "KW-1134", "KW-0998", "KW-0964"  # TM beta strand; Cell outer membrane; Secreted
# UniProt Domain names that are a whole membrane-spanning beta-barrel domain. How the list
# was chosen (no test-period data was used):
# - from the training period only (2023-07-01 to 2024-08-31): the Domain names on proteins
#   UniProt marks KW-1134 or KW-0998 there are "TBDR beta-barrel", "TonB-dependent
#   receptor-like beta-barrel", "Autotransporter" and "Porin";
# - from the PROSITE/Pfam/InterPro definitions of outer-membrane barrel domains, not from
#   any data: "YadA-like C-terminal membrane anchor" (IPR005594), the secretin barrel
#   ("secretion system secretin-like", IPR004846) and "Outer membrane protein
#   beta-barrel" (PF13505). (The secretin's N-terminal domains are periplasmic.)
# Every name is anchored to the whole domain name. Region notes such as aerolysin's "Part
# of the transmembrane beta-barrel after ... insertion into the host membrane" describe a
# conditional pore, not the protein's own barrel domain, and are not counted. Many soluble
# domains are barrels too (TIM barrel, RNA helicase aquarius beta-barrel), and the plug,
# POTRA and OmpA-like domains of outer-membrane proteins sit outside the membrane.
BARREL_DOMAIN = re.compile(r"^(?:(?:TBDR|TonB-dependent receptor(?:-like)?) beta-barrel|porin|autotransporter|"
                           r"YadA-like C-terminal membrane anchor|(?:type II/III )?secretion system secretin-like(?: domain)?|"
                           r"outer membrane protein beta-barrel(?: domain)?)$", re.I)
# A TonB-dependent receptor's barrel follows its plug. When UniProt names only the plug,
# the BARREL_AFTER_PLUG residues after it are taken as (the start of) the barrel, for
# proteins UniProt marks as outer-membrane or beta-strand (KW-0998, KW-1134). The
# training period has such plug-only proteins too ("TonB-dependent receptor plug").
TBDR_PLUG = re.compile(r"^(?:TBDR|TonB-dependent receptor) plug(?: domain)?$", re.I)
BARREL_AFTER_PLUG = 300
# The curators' rule: a structure is a membrane protein when a significant membrane-crossing
# part is resolved in it. Chosen on the training period only (see reports/membrane.md,
# "Resolved membrane-crossing segments"): structures where one protein UniProt marks
# transmembrane has at least CORE_MIN_SEGMENTS resolved segments of CORE_MIN_LENGTH or more
# residues, and that protein is at least CORE_MIN_SHARE of the residues modelled (a core
# protein, not a minor partner), were membrane proteins in 95.7% of the 1,823 training
# structures with three or more, and in 76% (31 of 41) of those the model scored below 50%
# (held-out predictions, counting barrels too). One or two segments were not enough (22% and
# 38% of those scored below 50%): pili, amyloid fibrils, phage coats, viral envelopes and
# helices UniProt marks in soluble enzymes. Segments are counted in one chain: a viral
# polyprotein's envelope proteins are separate chains with one or two anchors each. A
# resolved membrane beta-barrel (a barrel domain, or CORE_MIN_STRANDS resolved transmembrane
# strands of a beta-strand protein) counts too. Such structures are shown at least TM_FLOOR
# (from the 76% above; training period only), so family statistics or "secreted"/"cytoplasm"
# keywords cannot pull them below 50%. The rule is only a floor, not also a weight the model
# learns: as a weight it took strength from tm_modelled and lowered membrane complexes whose
# segments are spread over small subunits (signal peptidase, SecYEG with SecA, BAM).
CORE_MIN_SEGMENTS, CORE_MIN_LENGTH, CORE_MIN_SHARE, CORE_MIN_STRANDS = 3, 15, 0.25, 8
TM_FLOOR = 0.7
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


def _fetch(ids: list[str], chains: bool = False):
    """chains: take the coverage from rcsb.entry_chains (every chain and the sequences),
    to be checked with checked_coverage, instead of rcsb.entry_coverage. Used when scoring
    new structures; training reads the PDB's alignment as it is (see checked_coverage)."""
    ext, fam, det = rcsb.entry_extras(ids), rcsb.entry_families(ids), rcsb.entry_details(ids)
    cov = rcsb.entry_chains(ids) if chains else rcsb.entry_coverage(ids)
    # Segments for every linked protein, not only the aligned ones, so "UniProt has
    # no segments for it" is known rather than guessed.
    accs = sorted({a for v in cov.values() for e in (v or []) for a in e.get("align") or {}} |
                  {u for d in det.values() for u in _uniprots(d)})
    return ext, fam, det, cov, uniprot.membrane_features(accs)


def barrel_domains(named: list, keywords=()) -> list[list[int]]:
    """[[first, last], ...] of UniProt Domain/Region features that are a membrane-spanning
    beta-barrel (e.g. 'TBDR beta-barrel', 'Autotransporter'). Barrel-shaped domains that
    stay outside the membrane (OmpA-like, TIM barrel, plug, POTRA, ...) are not counted.

    keywords: the protein's UniProt keywords. If it is an outer-membrane or beta-strand
    protein with a TonB-dependent receptor plug but no named barrel, the residues after
    the plug are taken as the barrel."""
    out = [[b, e] for b, e, name in named or [] if BARREL_DOMAIN.search(name.strip())]
    if not out and {KW_OUTER, KW_STRAND} & set(keywords):
        out = [[e + 1, e + BARREL_AFTER_PLUG] for b, e, name in named or [] if TBDR_PLUG.search(name.strip())]
    return out


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


def _to_entity(u: int, regions: list) -> int | None:
    """UniProt residue u in the chain's numbering (None if it is not in the construct)."""
    for eb, rb, n in regions:
        if rb <= u < rb + n:
            return eb + u - rb
    return None


def tm_by_protein(cov: list | None, tmseg: dict[str, list]) -> tuple[dict[str, dict], int]:
    """({accession: {"segments": [[first, last], ...] resolved, "residues": modelled}}, total).

    A segment is resolved as in tm_coverage (in the construct, at least half of it
    modelled). Both are taken from one chain (entity) of the protein, never added up
    over chains: "segments" are those of its chain with the most resolved segments
    (of CORE_MIN_LENGTH or more, then of any length), "residues" the modelled residues
    of its largest chain. A viral polyprotein whose products are separate chains (prM
    and E, E1 and E2, each with one or two anchors) is one UniProt entry, but no chain
    of it crosses the membrane three times. total sums "residues" over the structure's
    UniProt-linked proteins, so the share of a protein tells a core protein from a
    minor partner."""
    out: dict[str, dict] = {}
    long_ = lambda segs: (sum(1 for b, e in segs if e - b + 1 >= CORE_MIN_LENGTH), len(segs))
    for ent in cov or []:
        unobs = ent.get("unobs") or []
        hidden = lambda q: any(lo <= q <= hi for lo, hi in unobs)
        for acc, regions in (ent.get("align") or {}).items():
            seen = sum(1 for eb, _, n in regions for q in range(eb, eb + n) if not hidden(q))
            d = out.setdefault(acc, {"segments": [], "residues": 0})
            d["residues"] = max(d["residues"], seen)
            segs = []
            for b, e in tmseg.get(acc) or []:
                pos = [q for q in (_to_entity(u, regions) for u in range(b, e + 1)) if q is not None]
                if len(pos) >= (e - b + 1) / 2 and sum(1 for q in pos if not hidden(q)) >= (e - b + 1) / 2 \
                        and [b, e] not in segs:
                    segs.append([b, e])
            if long_(segs) > long_(d["segments"]):
                d["segments"] = segs
    return out, sum(d["residues"] for d in out.values())


def core_tm(by_protein: dict[str, dict], total: int, tm_accessions, strand_accessions=()) -> dict:
    """The protein that carries the curators' rule, if any.

    {"protein": accession or None, "segments": its resolved segments of CORE_MIN_LENGTH or
    more, "share": its share of the modelled residues, "strands": resolved transmembrane
    strands of a beta-strand protein (any length; UniProt strands are 7-12 residues),
    "resolved": its resolved segments of any length (in the same chain),
    "partner_segments": the most such segments in a minor partner (under CORE_MIN_SHARE)}."""
    best = {"protein": None, "segments": 0, "resolved": 0, "share": 0.0, "strands": 0, "partner_segments": 0}
    for acc, d in by_protein.items():
        share = d["residues"] / total if total else 0.0
        if acc in strand_accessions:
            best["strands"] = max(best["strands"], len(d["segments"]) if share >= CORE_MIN_SHARE else 0)
        if acc not in tm_accessions:
            continue
        n = sum(1 for b, e in d["segments"] if e - b + 1 >= CORE_MIN_LENGTH)
        if share < CORE_MIN_SHARE:  # a minor partner (reported, but not the rule)
            best["partner_segments"] = max(best["partner_segments"], n)
            continue
        if (n, share) > (best["segments"], best["share"]):
            best.update(protein=acc, segments=n, resolved=len(d["segments"]), share=round(share, 3))
    return best


ALIGN_TRUST = 0.98   # identity over the PDB's aligned regions at or above which they are used as they are
ALIGN_MIN = 0.7      # ... below which (and if the sequences cannot be lined up either) the chain is not checked
ALIGN_BLOCK = 6      # shortest identical stretch that anchors a re-alignment
ALIGN_FILL = 40      # longest run of substitutions bridged between two anchors on the same diagonal
ALIGN_MAX_CELLS = 40_000_000  # entity length x UniProt length above which no re-alignment is tried


def _identical(regions: list, seq: str, ref: str) -> tuple[int, int]:
    """(aligned residues, of those identical in the two sequences)."""
    n = same = 0
    for eb, rb, k in regions:
        for j in range(k):
            n += 1
            q, u = eb + j, rb + j
            same += 1 <= q <= len(seq) and 1 <= u <= len(ref) and seq[q - 1] == ref[u - 1]
    return n, same


def _clip(regions: list, ref_len: int) -> list:
    """Regions without the residues numbered past the end of the UniProt sequence."""
    out = []
    for eb, rb, k in regions:
        k = min(k, ref_len - rb + 1)
        if k > 0 and rb >= 1:
            out.append([eb, rb, k])
    return out


@functools.lru_cache(maxsize=4096)
def _seq_blocks(seq: str, ref: str) -> tuple:
    return tuple((b.a, b.b, b.size) for b in difflib.SequenceMatcher(None, seq, ref, autojunk=False).get_matching_blocks()
                 if b.size >= ALIGN_BLOCK)


def seq_regions(seq: str, ref: str) -> list[list[int]]:
    """[[entity_start, uniprot_start, length], ...] from the sequences themselves:
    identical stretches of ALIGN_BLOCK or more, in order, joined across substitutions."""
    if not seq or not ref or len(seq) * len(ref) > ALIGN_MAX_CELLS:
        return []
    out: list[list[int]] = []
    for a, b, size in _seq_blocks(seq, ref):
        if out:
            eb, rb, k = out[-1]
            gap = a + 1 - (eb + k)
            if gap == b + 1 - (rb + k) and 0 <= gap <= ALIGN_FILL:
                out[-1][2] = k + gap + size
                continue
        out.append([a + 1, b + 1, size])
    return out


def _unclaimed(new: list, claimed: set, regions: list) -> list:
    """The stretches of the re-derived regions on chain residues no PDB region claims
    (claimed) and on UniProt residues the PDB's regions for this protein do not reach."""
    reached = {rb + j for _, rb, k in regions for j in range(k)}
    free = [(eb + j, rb + j) for eb, rb, k in new for j in range(k) if eb + j not in claimed and rb + j not in reached]
    out: list[list[int]] = []
    for q, u in free:
        if out and out[-1][0] + out[-1][2] == q and out[-1][1] + out[-1][2] == u:
            out[-1][2] += 1
        else:
            out.append([q, u, 1])
    return [r for r in out if r[2] >= ALIGN_BLOCK]


def check_alignment(regions: list, seq: str, ref: str, claimed: set | None = None) -> tuple[list, str]:
    """(regions to use, how they were got) for one chain against one UniProt sequence.

    The PDB often gives one gap-free region although the construct has a deletion or an
    insertion, or numbers the chain by another isoform than UniProt's entry: every residue
    after that point is then shifted, and a membrane segment looked up there is the wrong
    piece of chain. It also leaves parts out (the stretch after a deleted loop). So the
    regions are compared with the two sequences:
    "pdb": they agree (identity ALIGN_TRUST or more) and lining up the sequences finds no
           more, used as given. claimed: the chain holds several UniProt proteins (a
           fusion or chimera) and these are the chain residues the PDB gives to any of
           them. Its division is then kept, since two related proteins (HCN1 and HCN2)
           would both line up with the same stretch, and only stretches it gives to no
           protein are added;
    "sequence": lining up the sequences explains more residues (or as many, where the
           PDB's regions do not agree with the sequences): the regions are re-derived;
    "pdb-homolog": a related sequence (identity ALIGN_MIN or more) that the re-alignment
           does not improve on, used as given;
    "unchecked": neither works, [] is returned: nothing may be read off this chain."""
    regions = _clip(regions, len(ref))
    n, same = _identical(regions, seq, ref)
    trusted = bool(n) and same >= ALIGN_TRUST * n
    new = seq_regions(seq, ref)
    if trusted and claimed is not None:
        more = _unclaimed(new, claimed, regions)
        return (sorted(regions + more), "sequence") if more else (regions, "pdb")
    new_same = _identical(new, seq, ref)[1]
    if new_same > same - (not trusted) and new_same >= min(20, 0.8 * len(seq)):
        return new, "sequence"
    if trusted:
        return regions, "pdb"
    if n and same >= ALIGN_MIN * n:
        return regions, "pdb-homolog"
    return [], "unchecked"


def _observed(regions: list, unobs: list) -> int:
    return sum(1 for eb, _, n in regions for q in range(eb, eb + n) if not any(lo <= q <= hi for lo, hi in unobs))


def _resolved(seg: list, regions: list, unobs: list) -> bool:
    """A UniProt segment is resolved in a chain: in the construct and at least half modelled."""
    b, e = seg
    pos = [q for q in (_to_entity(u, regions) for u in range(b, e + 1)) if q is not None]
    return len(pos) >= (e - b + 1) / 2 and \
        sum(1 for q in pos if not any(lo <= q <= hi for lo, hi in unobs)) >= (e - b + 1) / 2


def checked_coverage(item: list | None, feats: dict | None = None, kw: dict | None = None,
                     proteins: dict | None = None) -> list | None:
    """An rcsb.entry_chains item as an rcsb.entry_coverage item that can be trusted.

    Per entity: {"align": regions per accession after check_alignment (an accession that
    could not be checked is left out), "unobs": unobserved residues of its best chain,
    "unchecked": [accessions left out], "realigned": [accessions whose regions were
    re-derived from the sequences]}.

    The best chain is the copy with the most resolved membrane segments or barrel domains
    (of CORE_MIN_LENGTH or more first), then the most modelled residues: copies of a chain
    are often modelled to different extents, and the first one may lack what another has.

    feats: uniprot.membrane_features output; kw: {accession: UniProt keywords};
    proteins: uniprot.entries output (without a UniProt sequence the PDB's regions stand)."""
    if item is None:
        return None
    feats, kw, proteins = feats or {}, kw or {}, proteins or {}
    out = []
    for ent in item:
        align, unchecked, realigned = {}, [], []
        given = ent.get("align") or {}
        claimed = {eb + j for regions in given.values() for eb, _, k in regions for j in range(k)} if len(given) > 1 else None
        for acc, regions in given.items():
            ref, seq = (proteins.get(acc) or {}).get("seq"), ent.get("seq")
            if not ref or not seq:
                if regions:
                    align[acc] = regions
                else:
                    unchecked.append(acc)
                continue
            regions, how = check_alignment(regions, seq, ref, claimed)
            if how == "unchecked":
                unchecked.append(acc)
            else:
                align[acc] = regions
                if how == "sequence":
                    realigned.append(acc)
        targets = {}
        for acc in align:
            f = feats.get(acc) or {}
            targets[acc] = f.get("tm") or barrel_domains(f.get("named") or [], kw.get(acc, ()))
        def worth(unobs):
            hits = [e - b + 1 >= CORE_MIN_LENGTH for acc, regions in align.items() for b, e in targets[acc]
                    if _resolved([b, e], regions, unobs)]
            return sum(hits), len(hits), sum(_observed(regions, unobs) for regions in align.values()), \
                -sum(hi - lo + 1 for lo, hi in unobs)
        chains = ent.get("chains") or [ent.get("unobs") or []]
        out.append({"align": align, "unobs": max(chains, key=worth), "unchecked": unchecked, "realigned": realigned})
    return out


TABLE_MAX_SEGMENTS = 30  # segmentRanges kept per protein (titin-sized lists are no help)


def _ranges(nums) -> list[list[int]]:
    """Sorted residue numbers as merged [[first, last], ...]."""
    out: list[list[int]] = []
    for n in sorted(set(nums)):
        if out and n <= out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return out


def _range_text(ranges: list, gap: int = 25, most: int = 4) -> str:
    """'168-628' or '31-200, 260-489': ranges for a sentence, short loops not shown as breaks."""
    merged: list[list[int]] = []
    for b, e in ranges:
        if merged and b - merged[-1][1] <= gap:
            merged[-1][1] = e
        else:
            merged.append([b, e])
    text = ", ".join(f"{b}-{e}" for b, e in merged[:most])
    return text + (f" and {len(merged) - most} more" if len(merged) > most else "")


def _chain_view(ent: dict, acc: str, segs: list) -> dict:
    """One chain (entity) of a protein: modelled UniProt residues, and which segments are
    resolved (in the construct and at least half modelled, as tm_coverage counts them)."""
    unobs, regions = ent.get("unobs") or [], (ent.get("align") or {}).get(acc) or []
    hidden = lambda q: any(lo <= q <= hi for lo, hi in unobs)
    modelled = _ranges(rb + q - eb for eb, rb, n in regions for q in range(eb, eb + n) if not hidden(q))
    return {"modelled": modelled, "flags": [_resolved(g, regions, unobs) for g in segs],
            "residues": sum(e - b + 1 for b, e in modelled), "realigned": acc in (ent.get("realigned") or [])}


def _within(ranges: list, b: int, e: int) -> list[list[int]]:
    return [[max(lo, b), min(hi, e)] for lo, hi in ranges if hi >= b and lo <= e]


def _barrel_name(named: list, rng: list) -> str:
    for b, e, name in named or []:
        if [b, e] == rng and BARREL_DOMAIN.search(name.strip()):
            return name.strip()
    return "beta-barrel after the TBDR plug (estimated)"


def tm_table(det: dict | None, cov: list | None, feats: dict | None, kw: dict | None = None,
             chains: dict | None = None, proteins: dict | None = None) -> list[dict]:
    """Per-protein evidence for one structure, so a curator can check the right protein.

    One row per UniProt protein (and one per protein entity with no UniProt link), from
    its chain with the most resolved segments, as tm_by_protein picks it:
    {"uniprot", "name" (UniProt's protein name; the PDB entity description if unknown),
     "chains" (count), "length" (UniProt sequence length, 0 if not known),
     "modelled": [[first, last], ...] observed residues in UniProt numbering,
     "segments": number of UniProt transmembrane/intramembrane segments,
     "resolved": how many of them are resolved (same criterion as the score),
     "segmentRanges": [[first, last, resolved], ...] (at most TABLE_MAX_SEGMENTS),
     "barrel": {"barrel": name, "range": [first, last], "resolved": bool} when UniProt maps
               no segments but a membrane-spanning barrel domain (barrel_domains),
     "status": "resolved" (all), "partial", "unresolved", "none" (no membrane part in
               UniProt) or "unknown" (could not be checked: no UniProt link, UniProt has
               deleted the entry ("gone": True), or the chain's sequence could not be lined
               up with UniProt's; resolved flags are then None and "modelled" is empty, so
               no residue numbers are shown that are not UniProt's),
     "pdbName": what the PDB calls the chain, when that is not UniProt's name (a viral
               polyprotein is one UniProt entry; the PDB names the product, e.g. 'M protein'),
     "fusion": True if the PDB chain it is in also holds another UniProt protein,
     "realigned": True if the positions come from the sequences because the PDB's own
               residue mapping was shifted or incomplete (check_alignment),
     "note": the same in plain English}.
    Rows: proteins with resolved segments first (most first), then with unresolved ones,
    then without any, then those that could not be checked.

    det: rcsb.entry_details item; cov: checked_coverage output (or an rcsb.entry_coverage
    item); feats: uniprot.membrane_features output; kw: {accession: UniProt keywords};
    chains: {accession or entity description: number of chains}; proteins: uniprot.entries
    output ({accession: {"name", "seq"}}, None for an entry UniProt no longer has)."""
    feats, kw, chains, proteins = feats or {}, kw or {}, chains or {}, proteins or {}
    names: dict[str, str] = {}
    fused: set[str] = set()
    unlinked: list[str] = []
    for e in (det or {}).get("polymer_entities") or []:
        if (e.get("entity_poly") or {}).get("rcsb_entity_polymer_type") not in (None, "Protein"):
            continue
        desc = ((e.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()
        ups = (e.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []
        for u in ups:
            names.setdefault(u, desc)
        if len(ups) > 1:
            fused.update(ups)
        if not ups and desc not in unlinked:
            unlinked.append(desc)
    order = list(dict.fromkeys([a for ent in cov or [] for a in list(ent.get("align") or {}) + list(ent.get("unchecked") or [])]
                               + list(names)))
    rows = []
    for acc in order:
        f = feats.get(acc) or {}
        known = proteins.get(acc)
        gone = acc in proteins and known is None
        segs = sorted([list(s) for s in {tuple(s) for s in f.get("tm") or []}])
        barrels = [] if segs else barrel_domains(f.get("named") or [], kw.get(acc, ()))
        targets = segs or barrels
        long_ = lambda v: (sum(1 for (b, e), ok in zip(targets, v["flags"]) if ok and e - b + 1 >= CORE_MIN_LENGTH),
                           sum(v["flags"]), v["residues"])
        views = [] if gone else [_chain_view(ent, acc, targets) for ent in cov or [] if (ent.get("align") or {}).get(acc)]
        best = max(views, key=long_) if views else None
        row = {"uniprot": acc, "name": (known or {}).get("name") or names.get(acc, ""), "chains": chains.get(acc) or 1,
               "length": len((known or {}).get("seq") or ""), "modelled": best["modelled"] if best else [],
               "segments": len(segs), "resolved": 0, "segmentRanges": []}
        if acc in fused:
            row["fusion"] = True
        flags = best["flags"] if best else [None] * len(targets)
        where = _range_text(row["modelled"])
        only = f"only residues {where} were modelled" if where else "none of it was modelled"
        what = (f"{len(segs)} membrane-crossing segment{'' if len(segs) == 1 else 's'} ({_range_text(segs, gap=0, most=3)})"
                if segs else f"a beta-barrel domain ({_range_text(barrels, gap=0)})" if barrels else "")
        if segs:
            n, k = len(segs), sum(map(bool, flags))
            row["resolved"] = k
            row["segmentRanges"] = [[b, e, ok] for (b, e), ok in zip(segs, flags)][:TABLE_MAX_SEGMENTS]
            missing = [f"{b}-{e}" for (b, e), ok in zip(segs, flags) if not ok]
            miss = ", ".join(missing[:3]) + (f" and {len(missing) - 3} more" if len(missing) > 3 else "")
            if k == n:
                row["status"] = "resolved"
                row["note"] = (f"its 1 membrane-crossing segment ({segs[0][0]}-{segs[0][1]}) is in the structure" if n == 1
                               else f"{'both' if n == 2 else 'all ' + str(n)} membrane-crossing segments are in the structure")
            elif k:
                row["status"] = "partial"
                row["note"] = f"{k} of its {n} membrane-crossing segments {'is' if k == 1 else 'are'} in the structure (not {miss})"
            else:
                row["status"] = "unresolved"
                row["note"] = (f"its 1 membrane-crossing segment ({miss}) is not in the structure: {only}" if n == 1
                               else f"none of its {n} membrane-crossing segments ({miss}) is in the structure: {only}")
        elif barrels:
            i = max(range(len(barrels)), key=lambda j: bool(flags[j]))
            b, e = barrels[i]
            row["barrel"] = {"barrel": _barrel_name(f.get("named") or [], barrels[i]), "range": [b, e], "resolved": flags[i]}
            row["status"] = "resolved" if flags[i] else "unresolved"
            inside = _within(row["modelled"], b, e)
            share = sum(hi - lo + 1 for lo, hi in inside) / (e - b + 1)
            row["note"] = (f"beta-barrel domain {b}-{e} is not in the structure: {only}" if not flags[i]
                           else f"beta-barrel domain {b}-{e} is in the structure" if share >= 0.9
                           else f"{'most' if share >= 0.75 else 'part'} of beta-barrel domain {b}-{e} is in the structure "
                                f"(modelled {_range_text(inside)})")
        else:
            row["status"] = "none"
            row["note"] = "no membrane-crossing segments in UniProt"
        if gone:
            # UniProt answers nothing for a deleted entry, which is not "it has no segments".
            row.update(status="unknown", gone=True,
                       note="UniProt no longer has this entry (it was deleted or merged into another), "
                            "so its membrane segments could not be checked")
        elif not best and targets:
            # No residue numbers are shown for this chain: they would not be UniProt's.
            row["status"] = "unknown"
            row["note"] = (f"UniProt maps {what}, but this chain's sequence could not be lined up with the UniProt entry, "
                           f"so {'it' if len(targets) == 1 else 'they'} could not be checked automatically")
        elif best and best["realigned"]:
            row["realigned"] = True
            if targets:
                row["note"] += " (positions lined up by sequence: the PDB's own residue mapping to UniProt is shifted or incomplete for this chain)"
        if acc not in fused and known and names.get(acc) and names[acc].lower() != row["name"].lower():
            row["pdbName"] = names[acc] if len(names[acc]) <= 80 else names[acc][:77].rstrip() + "..."
        rows.append(row)
    for desc in unlinked:
        rows.append({"uniprot": None, "name": desc, "chains": chains.get(desc) or 1, "length": 0, "modelled": [],
                     "segments": 0, "resolved": 0, "segmentRanges": [], "status": "unknown",
                     "note": "not linked to UniProt, so its membrane segments could not be checked"})
    rank = {"resolved": 0, "partial": 0, "unresolved": 1, "none": 2, "unknown": 3}
    done = lambda r: r["resolved"] or int(bool((r.get("barrel") or {}).get("resolved")))
    rows.sort(key=lambda r: (rank[r["status"]], -done(r)))
    return rows


def chain_counts(ent: dict | None) -> dict[str, int]:
    """{accession (or, for entities with no UniProt link, entity description): chains} from
    an rcsb.entry_entities item."""
    out: dict[str, int] = defaultdict(int)
    for e in (ent or {}).get("polymer_entities") or []:
        ids = e.get("rcsb_polymer_entity_container_identifiers") or {}
        n = len(ids.get("auth_asym_ids") or []) or 1
        for key in ids.get("uniprot_ids") or [((e.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()]:
            out[key] += n
    return dict(out)


def name_proteins(table: list[dict], most: int = 3) -> str:
    """'A (P1) 6, B (P2) 4 and 1 more': the proteins whose membrane-crossing segments are
    resolved, with how many each. Words that all listed names start with are said once
    (not when a name has a comma, as the PDB writes fusion chains: 'A,B chimera')."""
    rows = [r for r in table if r.get("uniprot") and r.get("resolved")]
    shown = rows[:most]
    words = [(r.get("name") or "").split() for r in shown]
    common = 0
    while len(shown) > 1 and not any("," in (r.get("name") or "") for r in shown) and all(len(w) > common + 1 for w in words) and len({w[common].lower() for w in words}) == 1:
        common += 1
    parts = []
    for i, (r, w) in enumerate(zip(shown, words)):
        name = " ".join(w if i == 0 else w[common:])
        name = name if len(name) <= 48 else name[:45].rstrip() + "..."
        parts.append((f"{name} ({r['uniprot']})" if name else r["uniprot"]) + f" {r['resolved']}")
    return ", ".join(parts) + (f" and {len(rows) - most} more" if len(rows) > most else "")


def raw_facts(pdb: str, ext: dict | None, fam: dict | None, det: dict | None,
              cov: list | None = None, feats: dict | None = None) -> dict | None:
    """Release-day facts about a structure (no labels, no PDB annotations).

    feats: uniprot.membrane_features output ({acc: {"tm": [...], "named": [...]}}).
    """
    if not det:
        return None
    tmseg = {a: v["tm"] for a, v in feats.items()} if feats is not None else None
    kw = (ext or {}).get("kw") or {}
    # Named barrel domains count whatever the keywords: unreviewed porins often carry no
    # transmembrane keyword (8XWO chitoporin), and the names are anchored to whole domains.
    barrel = ({a: barrel_domains(v["named"], kw.get(a, ())) for a, v in feats.items()}
              if feats is not None else None)
    prot = [e for e in det.get("polymer_entities") or []
            if (e.get("entity_poly") or {}).get("rcsb_entity_polymer_type") == "Protein"]
    ups = _uniprots(det) - FUSION_UNIPROT
    tm_ups = {u for u in ups if TM_KEYWORDS & set(kw.get(u, ()))}
    all_kw = {k for u in ups for k in kw.get(u, ())}
    ligs = {((n.get("nonpolymer_comp") or {}).get("chem_comp") or {}).get("id") for n in det.get("nonpolymer_entities") or []}
    descs = [((e.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "") for e in prot]
    ipr = (fam or {}).get("ipr") or {}
    has_tm_protein = any(tmseg.get(u) for u in ups) if tmseg is not None else None
    has_barrel = any(barrel.get(u) for u in ups) if barrel is not None else None
    # Outer-membrane or beta-strand proteins UniProt calls transmembrane but maps neither
    # segments nor a barrel domain for, so whether their membrane part is resolved cannot
    # be checked. Secreted proteins (toxins UniProt gives the transmembrane keyword, such
    # as botulinum NTNH or latrotoxin) are not counted: their pore forms only in a host.
    tm_unmapped = (sorted(u for u in tm_ups if not tmseg.get(u) and not barrel.get(u)
                          and {KW_OUTER, KW_STRAND} & set(kw.get(u, ())) and KW_SECRETED not in kw.get(u, ()))
                   if tmseg is not None else None)
    # Families of the proteins UniProt places in a membrane, so a soluble partner's family
    # (a ubiquitin-like or TPR protein bound to a barrel) does not speak for the complex.
    mem_ups = {u for u in ups if (TM_KEYWORDS | {KW_OUTER, KW_STRAND}) & set(kw.get(u, ()))}
    tm_total, tm_modelled = tm_coverage(cov, tmseg or {}) if cov is not None else (None, None)
    # Membrane-spanning beta-barrel domains, mapped onto the structure the same way.
    bar_total, bar_modelled = tm_coverage(cov, barrel or {}) if cov is not None else (None, None)
    core = None
    if cov is not None and tmseg is not None:
        strand_ups = {u for u in ups if KW_STRAND in kw.get(u, ())}
        core = core_tm(*tm_by_protein(cov, tmseg), tm_ups, strand_ups)
    return {
        "core_tm": core,
        "has_tm_protein": has_tm_protein, "has_barrel": has_barrel, "tm_unmapped": tm_unmapped,
        "tm_in_construct": tm_total, "tm_modelled": tm_modelled,
        "barrel_in_construct": bar_total, "barrel_modelled": bar_modelled,
        "pdb": pdb, "released": (ext or {}).get("released", ""),
        "uniprots": sorted(ups), "tm_uniprots": sorted(tm_ups), "keywords": sorted(all_kw),
        "detergents": sorted(ligs & DETERGENTS), "lipids": sorted(ligs & LIPIDS),
        "nanodisc": any(NANODISC.search(d) for d in descs),
        "ipr": sorted({x for u in ups for x in ipr.get(u, ())}),
        "ipr_membrane": sorted({x for u in mem_ups for x in ipr.get(u, ())}),
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
    # A membrane protein here means one whose membrane-spanning part is in the structure:
    # resolved helical segments (UniProt transmembrane features) or a resolved
    # membrane beta-barrel (UniProt transmembrane strands or a barrel domain).
    if f.get("tm_modelled") is not None:
        bar = f.get("barrel_modelled") or 0
        x["tm_modelled"] = float(f["tm_modelled"] > 0 or bar > 0)
        x["barrel_modelled"] = float(bar > 0)
        x["tm_missing"] = float((bool(f.get("has_tm_protein")) or bool(f.get("has_barrel")))
                                and f["tm_modelled"] == 0 and bar == 0)
    # The curators' rule (see TM_FLOOR): not a weight the model learns but a floor on the
    # probability, so keys starting "rule:" are left out of fitting.
    if f.get("core_tm") is not None:
        core, bar = f["core_tm"], (f.get("barrel_modelled") or 0) > 0
        strands = core.get("strands", 0) >= CORE_MIN_STRANDS
        x["rule:tm_resolved"] = float(core["segments"] >= CORE_MIN_SEGMENTS or strands or bar)
        why["rule"] = {"protein": core["protein"], "segments": core["segments"],
                       "resolved": core.get("resolved", core["segments"]), "barrel": bar or strands}
    # UniProt calls it transmembrane but maps nothing to check: unknown, not "missing".
    if f.get("tm_unmapped") is not None:
        x["tm_unmapped"] = float(bool(f["tm_unmapped"]) and not f.get("has_tm_protein") and not f.get("has_barrel"))
    x["kw_strand"] = float(KW_STRAND in kw)
    x["kw_outer"] = float(KW_OUTER in kw)
    tcs = sorted({tc for u in f["uniprots"] for tc in ctx.tcdb.get(u, ())})
    x["tcdb_transporter"] = float(any(tc[0] in "123" for tc in tcs))
    x["tcdb_other"] = float(bool(tcs) and not x["tcdb_transporter"])
    x["tcdb_barrel"] = float(any(tc.startswith("1.B.") for tc in tcs))  # outer-membrane beta-barrel porins
    why["tcdb"] = sorted(tcs, key=lambda tc: not tc.startswith("1.B."))[:2] if x["tcdb_barrel"] else tcs[:2]
    x["detergent"] = float(bool(f["detergents"]))
    x["lipid"] = float(bool(f["lipids"]))
    x["nanodisc"] = float(f["nanodisc"])
    earlier = [u for u in f["uniprots"] if any(d and d < f["released"] for d in ctx.protein_dates.get(u, ()))]
    x["same_protein"] = float(bool(earlier))
    why["same_protein"] = ctx.protein_example.get(earlier[0]) if earlier else None
    best, best_ipr = 0.0, None
    # Family of a membrane protein in it if there is one, else of any protein in it.
    for i in f.get("ipr_membrane") or f["ipr"]:
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
    keys = sorted(k for k in counts if not k.startswith("rule:") and (not k.startswith("t:") or counts[k] >= min_count))
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
    """Calibrated probability, raised to model["tm_floor"] when the curators' rule holds."""
    p = predict_model(x, model)
    floor = model.get("tm_floor")
    return max(p, floor) if floor and x.get("rule:tm_resolved") else p


def predict_model(x: dict[str, float], model: dict) -> float:
    """The calibrated model alone, without the curators' rule."""
    z = model["bias"] + sum(model["weights"].get(k, 0.0) * v for k, v in x.items())
    p = 1 / (1 + math.exp(-z))
    if model.get("isotonic"):  # observed frequencies on held-out training structures
        return _apply_isotonic(model["isotonic"], p)
    cal = model.get("calibration")
    if cal:  # Platt scaling fitted on held-out weeks
        z2 = cal["a"] * math.log(max(p, 1e-6) / max(1 - p, 1e-6)) + cal["b"]
        p = 1 / (1 + math.exp(-z2))
    return p


def explain(x: dict[str, float], why: dict, f: dict, p: float, model: dict, n: int = 3,
            table: list[dict] | None = None) -> list[str]:
    """Plain-English reasons, strongest first.

    High score: the evidence for it, then "but" and the strongest evidence
    against. When the curators' rule holds (a significant membrane-crossing part
    of a main protein is resolved), that is always the first reason. Low score: the evidence against it, with the strongest evidence
    for it folded into the last reason as "although ..." (or, if nothing
    counts against it, "the evidence for it is limited: ..."). So "but" only
    ever introduces evidence against a membrane protein.

    table: tm_table rows. With more than one protein in the structure, the resolved
    segments are named by protein, so a curator opens the right UniProt entry. A reason
    that lists proteins ends in a count ("... (P0DH92) 4 and 3 more"); the site shows it
    as a sentence of its own.
    """
    w = model["weights"]
    reasons = []  # (signed contribution, text)
    rule_reason = None  # under the curators' rule the resolved part is the first reason
    def add(key, val, text):
        if text:
            reasons.append((w.get(key, 0.0) * val, text))
    if x.get("kw_tm"):
        add("kw_tm", 1, "UniProt lists " + ("all its proteins" if x.get("kw_tm_all") and len(f["uniprots"]) > 1
                                            else "a protein in it") + " as crossing the membrane")
    elif f["uniprots"]:
        reasons.append((-abs(w.get("kw_tm", 1.0)), "UniProt does not list any of its proteins as crossing the membrane"))
    if x.get("kw_strand"):
        add("kw_strand", 1, "UniProt says a protein in it crosses the membrane as beta strands, as in a beta-barrel")
    if x.get("kw_outer"):
        add("kw_outer", 1, "UniProt places a protein in it in the outer membrane of a bacterium")
    if x.get("kw_membrane") and not x.get("kw_tm"):
        add("kw_membrane", 1, "UniProt places a protein in it in a membrane")
    for k, place in KW_SOLUBLE.items():
        if x.get(f"kw_{k}"):
            add(f"kw_{k}", 1, f"UniProt says a protein in it is {place}")
    if x.get("tcdb_barrel"):  # one sentence, carrying both TCDB weights (1.B is also a transporter class)
        reasons.append((w.get("tcdb_barrel", 0.0) + w.get("tcdb_transporter", 0.0) * x.get("tcdb_transporter", 0.0),
                        f"TCDB classifies a protein in it as an outer-membrane barrel protein (TC {why['tcdb'][0]})"))
    elif x.get("tcdb_transporter"):
        add("tcdb_transporter", 1, f"TCDB classifies a protein in it as a channel or transporter (TC {why['tcdb'][0]})")
    elif x.get("tcdb_other"):
        add("tcdb_other", 1, f"TCDB lists a protein in it (TC {why['tcdb'][0]})")
    if x.get("tm_modelled"):
        parts = []
        if f.get("tm_modelled"):
            k = f["tm_modelled"]
            parts.append(f"{k} membrane-crossing {'segment is' if k == 1 else 'segments are'} resolved in the structure")
        if x.get("barrel_modelled"):
            parts.append("the beta-barrel that crosses the membrane is resolved in the structure")
        rule = why.get("rule") or {}
        named = name_proteins(table) if table and sum(1 for r in table if r.get("uniprot")) > 1 else ""
        if named and parts and f.get("tm_modelled"):
            # The number said is the sum of the per-protein counts listed after it (each
            # protein's best chain), not tm_modelled, which counts every chain entity and
            # so counts a protein twice when the PDB splits it over two entities.
            holders = [r for r in table if r.get("uniprot") and r.get("resolved")]
            k = sum(r["resolved"] for r in holders)
            parts[0] = (f"{k} membrane-crossing {'segment is' if k == 1 else 'segments are'} resolved in the structure" +
                        (f", {'' if k == 1 else 'all '}in {named.rsplit(' ', 1)[0]}" if len(holders) == 1 else f": {named}"))
        elif x.get("rule:tm_resolved") and rule.get("segments", 0) >= CORE_MIN_SEGMENTS and parts and f.get("tm_modelled") \
                and len(f["uniprots"]) > 1:  # with one UniProt protein there is nothing to point to
            parts[0] += (f" ({rule['resolved']} of them in UniProt {rule['protein']}, one of its main proteins)"
                         if f["tm_modelled"] > rule["resolved"] else f" (in UniProt {rule['protein']}, one of its main proteins)")
        add("tm_modelled", 1, " and ".join(parts))
        if x.get("barrel_modelled"):
            reasons[-1] = (reasons[-1][0] + w.get("barrel_modelled", 0.0), reasons[-1][1])
        if x.get("rule:tm_resolved"):
            rule_reason = reasons[-1][1]
    if x.get("tm_missing"):
        add("tm_missing", 1, "the part that crosses the membrane is not resolved, so only an outside part of the protein is in the structure")
    if x.get("tm_unmapped"):
        add("tm_unmapped", 1, "UniProt has not marked where it crosses the membrane, so that part could not be checked")
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
    supporting = sorted([r for r in reasons if r[0] > 0.05], key=lambda r: -r[0])
    against = sorted([r for r in reasons if r[0] < -0.05], key=lambda r: r[0])
    if rule_reason and p >= 0.5:
        supporting.sort(key=lambda r: r[1] != rule_reason)
    if p >= 0.5:
        out = [t for _, t in supporting[:n]]
        if against:
            out.append("but " + against[0][1])
        return out
    out = [t for _, t in against[:n]]
    if supporting:
        if out:
            out[-1] += ", although " + supporting[0][1]
        else:
            out = ["the evidence for it is limited: " + " and ".join(t for _, t in supporting[:2])]
    return out


def _labels(start: str, end: str, mpstruc_ids: set[str], monotopic_ids: set[str] = frozenset(),
            monotopic: str = "positive") -> tuple[list[str], set[str], list[str]]:
    """(IDs used, membrane-protein set, all IDs released in the window).

    Left out: OPM-only monotopic/peripheral and peptide entries. monotopic says
    what to do with structures mpstruc lists only as monotopic (no transmembrane
    source lists them): "positive" counts them as membrane proteins, "exclude"
    leaves them out, "negative" counts them as not membrane proteins (no
    membrane-spanning part)."""
    ids = released_between(start, end)
    opm = membranedb.opm()
    tm = annotated(start, end, ["PDBTM", "MemProtMD"]) | (set(ids) & (mpstruc_ids - monotopic_ids)) | \
        {p for p in ids if opm.get(p) == "transmembrane"}
    mono = (set(ids) & monotopic_ids) - tm
    pos = tm | mono if monotopic == "positive" else tm
    unsure = {p for p in ids if opm.get(p) in ("monotopic/peripheral", "peptide")} - pos
    if monotopic == "exclude":
        unsure |= mono
    elif monotopic == "negative":
        unsure -= mono
    return [p for p in ids if p not in unsure], pos, ids


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


def mpstruc_groups(db) -> dict[str, str]:
    """{PDB ID: mpstruc group name} for entries and their related structures."""
    out: dict[str, str] = {}
    for g in db.groups:
        for sg in g.subgroups:
            for e in sg.parents:
                for x in [e, *e.daughters]:
                    out.setdefault(x.pdbCode.upper(), g.name)
                    for r in x.related:
                        out.setdefault(r.pdb.upper(), g.name)
    return out


MONOTOPIC = "MONOTOPIC MEMBRANE PROTEINS"
BARREL_GROUP = "TRANSMEMBRANE PROTEINS: BETA-BARREL"
OUTER = "outer membrane (UniProt)"
# Recent structures curators asked about; their scores are reported from the same run.
SPOTLIGHT = ("38FI", "38AQ", "9PTT", "35YU")
# Structures with no resolved membrane-crossing part, which should stay low under the rule.
SPOTLIGHT_LOW = ("9OS3", "30AZ", "37UU", "9ZDF", "9X9Y", "31OL", "33WM")


def false_positive_split(rows: list[tuple[str, float, int, dict]]) -> dict:
    """Test structures labelled "not a membrane protein" but scored 50% or more, split into
    those UniProt calls an outer-membrane beta-barrel protein (KW-0998 and KW-1134), which
    the label sources probably missed, and the rest."""
    fp = [(pdb, p, f) for pdb, p, y, f in rows if p >= 0.5 and not y]
    barrel = sorted(pdb for pdb, _, f in fp if {KW_OUTER, KW_STRAND} <= set(f["keywords"]))
    return {"total": len(fp), "outer_membrane_barrel": len(barrel), "other": len(fp) - len(barrel),
            "outer_membrane_barrel_ids": barrel}


def group_recall(preds: list[tuple[str, float, int | None, dict]], groups: dict[str, str]) -> dict:
    """Share scored at least 50%, by kind: every structure in each mpstruc group
    (whatever its training label), and membrane proteins UniProt places in the
    outer membrane (KW-0998). y is None for structures left out of the labels."""
    kinds: dict[str, list[float]] = defaultdict(list)
    for pdb, p, y, f in preds:
        if pdb in groups:
            kinds[groups[pdb]].append(p)
        if y and KW_OUTER in f["keywords"]:
            kinds[OUTER].append(p)
    return {k: {"n": len(v), "recall_pct": round(100 * sum(q >= 0.5 for q in v) / len(v), 1)}
            for k, v in sorted(kinds.items())}


def segment_bucket(f: dict) -> str | None:
    """Which row of the resolved-segment table a structure falls in (None: nothing resolved)."""
    core = f.get("core_tm")
    if not core or not (f.get("tm_modelled") or f.get("barrel_modelled")):
        return None
    n = lambda k: "3+" if k >= 3 else str(k)
    if core["segments"]:
        return f"main protein, {n(core['segments'])} segment{'s' if core['segments'] > 1 else ''}"
    if core["partner_segments"]:
        return f"minor partner only, {n(core['partner_segments'])} segment{'s' if core['partner_segments'] > 1 else ''}"
    if f.get("barrel_modelled") or core.get("strands", 0) >= CORE_MIN_STRANDS:
        return "beta-barrel (domain or strands)"
    return f"only segments shorter than {CORE_MIN_LENGTH} residues"


def rule_effect(rows: list[tuple[str, float, float, int, dict]]) -> dict:
    """How the curators' rule changes scores. rows: (pdb, model alone, with the rule, label, facts).

    "segments": share of structures labelled membrane proteins, by what is resolved;
    "lifted": structures the rule takes from below 50% to TM_FLOOR."""
    seg: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for _, _, _, y, f in rows:
        b = segment_bucket(f)
        if b:
            seg[b][0] += y
            seg[b][1] += 1
    applies = [(pdb, q, p, y) for pdb, q, p, y, f in rows if p > q]
    lifted = sorted((pdb, y) for pdb, q, p, y in applies if q < 0.5)
    return {"segments": {k: {"n": v[1], "membrane_pct": round(100 * v[0] / v[1], 1)} for k, v in sorted(seg.items())},
            "raised": len(applies), "lifted": len(lifted),
            "lifted_membrane_pct": round(100 * sum(y for _, y in lifted) / len(lifted), 1) if lifted else None,
            "lifted_ids": [pdb for pdb, _ in lifted], "lifted_not_labelled": [pdb for pdb, y in lifted if not y]}


def _family_stats(rows) -> dict[str, list[int]]:
    stats: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for f, y in rows:
        for i in set(f["ipr"]):
            stats[i][0 if y else 1] += 1
    return dict(stats)


def train_and_test(folds: int = 5, monotopic: str | None = None, write: bool = True,
                   tradeoff: dict | None = None) -> dict:
    """Train on TRAIN, calibrate on held-out folds of it, report on TEST.

    Family statistics are computed without the structure being scored (from
    the other time blocks of the training period), and the calibration maps
    out-of-fold predictions to observed frequencies, so the percentage shown
    is what was actually seen on structures the model had not learnt from.
    monotopic: how mpstruc's monotopic-only structures are labelled (see _labels).
    """
    monotopic = monotopic or MONOTOPIC_LABEL
    db = load(source_xml())
    ctx = Context(db)
    groups = mpstruc_groups(db)
    mono_ids = {p for p, g in groups.items() if g == MONOTOPIC}
    ids_tr, pos_tr, _ = _labels(*TRAIN, ctx.mpstruc_ids, mono_ids, monotopic)
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
    oof, oof_rows = [], []
    for k in range(folds):
        m = fit([r[1] for r in rows if r[5] != k], [r[3] for r in rows if r[5] != k])
        held = [r for r in rows if r[5] == k]
        oof += [(predict_model(r[1], m), r[3]) for r in held]
        oof_rows += held
    model = fit([r[1] for r in rows], [r[3] for r in rows])
    model["isotonic"] = _isotonic(oof)
    model["tm_floor"] = TM_FLOOR
    model["prior"] = prior
    model["fam_stats"] = {k: v for k, v in _family_stats(facts).items() if sum(v) >= 2}
    ctx.fam_stats = model["fam_stats"]

    model["monotopic_label"] = monotopic

    ids_te, pos_te, all_te = _labels(*TEST, ctx.mpstruc_ids, mono_ids, monotopic)
    used = set(ids_te)
    ext, fam, det, cov, seg = _fetch(all_te)
    rows_te, kinds = [], []
    for pdb in all_te:
        f = raw_facts(pdb, ext.get(pdb), fam.get(pdb), det.get(pdb), cov.get(pdb), seg)
        if f:
            x, why = features(f, ctx, prior)
            y = int(pdb in pos_te) if pdb in used else None
            kinds.append((pdb, predict(x, model), y, f))
            if pdb in used:
                rows_te.append((pdb, x, y, f, why))
    oof_cal = [_apply_isotonic(model["isotonic"], p) for p, _ in oof]
    oof_eval = [(max(q, TM_FLOOR) if r[1].get("rule:tm_resolved") else q, r[3]) for q, r in zip(oof_cal, oof_rows)]
    rule_tr = [(r[0], q, p, r[3], r[4]) for q, (p, _), r in zip(oof_cal, oof_eval, oof_rows)]
    rule_te = [(pdb, predict_model(x, model), predict(x, model), y, f) for pdb, x, y, f, _ in rows_te]
    misses = sorted(({"pdb": pdb, "p": round(p, 3), "title": f["title"][:90]} for pdb, p, _, f in kinds
                     if groups.get(pdb) == BARREL_GROUP and p < 0.5), key=lambda m: -m["p"])
    result = {"monotopic_label": monotopic, "train": {"window": TRAIN, **evaluate_preds(oof_eval)},
              "test": {"window": TEST, **evaluate(rows_te, model), "by_kind": group_recall(kinds, groups),
                       "false_positives": false_positive_split([(pdb, p, y, f) for pdb, p, y, f in kinds if y is not None]),
                       "beta_barrel_misses": misses,
                       "model_alone_at_50": evaluate_preds([(q, y) for _, q, _, y, _ in rule_te])["at_50"]},
              "tm_rule": {"floor": TM_FLOOR, "train": rule_effect(rule_tr), "test": rule_effect(rule_te)}}
    if not write:
        return result
    if tradeoff is not None:
        result["monotopic_tradeoff"] = {**tradeoff, monotopic: _brief(result)}
    MODEL_PATH.write_text(json.dumps(model))
    later = sorted(set(released_between("2025-07-01", "2025-12-31")) & ctx.mpstruc_ids)
    sc = score_structures(later, db, model)
    ps = [v["p"] for v in sc.values()]
    result["spotlight"] = score_structures(list(SPOTLIGHT + SPOTLIGHT_LOW), db, model)
    result["mpstruc_later"] = {"window": ["2025-07-01", "2025-12-31"], "n": len(ps),
                               "at_least_50_pct": round(100 * sum(p >= 0.5 for p in ps) / len(ps), 1) if ps else None,
                               "at_least_10_pct": round(100 * sum(p >= 0.1 for p in ps) / len(ps), 1) if ps else None}
    examples = sorted(rows_te, key=lambda r: -predict(r[1], model))
    result["examples"] = [{"pdb": r[0], "membrane": r[2], "p": round(predict(r[1], model), 3),
                           "why": explain(r[1], r[4], r[3], predict(r[1], model), model)} for r in examples[::max(1, len(examples) // 40)]]
    (ROOT / "reports" / "membrane.json").write_text(json.dumps(result, indent=1))
    (ROOT / "reports" / "membrane.md").write_text(report_md(result, model))
    return result


def _brief(result: dict) -> dict:
    t = result["test"]
    return {"n": t["n"], "membrane": t["membrane"], "auc": t["auc"], "at_50": t["at_50"], "by_kind": t["by_kind"]}


def load_model() -> dict | None:
    return json.loads(MODEL_PATH.read_text()) if MODEL_PATH.exists() else None


def _scored(ids: list[str], db=None, model: dict | None = None):
    """Yields (pdb, probability, reasons, tm_table rows) for each structure the PDB knows."""
    ctx = Context(db or load(source_xml()), model.get("fam_stats"))
    ext, fam, det, chains, seg = _fetch(ids, chains=True)
    ents = rcsb.entry_entities(ids)  # chains per protein
    try:
        proteins = uniprot.entries(sorted(seg))
    except Exception as err:  # then the PDB's alignment stands unchecked, and rows are named by the PDB
        print(f"membrane: UniProt names and sequences unavailable ({err})")
        proteins = {}
    for pdb in ids:
        # The same checked coverage for the score and for the table.
        cov = {pdb: checked_coverage(chains.get(pdb), seg, (ext.get(pdb) or {}).get("kw"), proteins)}
        f = raw_facts(pdb, ext.get(pdb), fam.get(pdb), det.get(pdb), cov.get(pdb), seg)
        if not f:
            continue
        x, why = features(f, ctx, model["prior"])
        p = predict(x, model)
        table = tm_table(det.get(pdb), cov.get(pdb), seg, (ext.get(pdb) or {}).get("kw"),
                         chain_counts(ents.get(pdb)), proteins)
        yield pdb, p, explain(x, why, f, p, model, table=table), table


def score_structures(ids: list[str], db=None, model: dict | None = None) -> dict[str, dict]:
    """{pdb: {"p": probability, "why": [reasons], "tmTable": tm_table rows}} for new structures."""
    model = model or load_model()
    if not model:
        return {}
    return {pdb: {"p": round(p, 3), "why": why, "tmTable": table} for pdb, p, why, table in _scored(ids, db, model)}


def backfill_table(path: Path | None = None, today: datetime.date | None = None) -> dict:
    """Add tmEvidence to every candidate in data/candidates.json without a weekly run.

    Reads the file, computes each candidate's per-protein table from the caches (and the
    network for anything missing), replaces membraneWhy only where its text changed and
    the score is still the one stored, and rewrites the file as automation.main writes it
    (compact separators), leaving everything else as it was, except "tmChecked": the day
    the tables were computed."""
    from .automation import SITE_DATA
    path = path or SITE_DATA / "candidates.json"
    text = path.read_text()
    data = json.loads(text)
    if json.dumps(data, separators=(",", ":")) != text.strip():
        print("membrane table: note: the file was not in automation.py's compact format; it will be now")
    model = load_model()
    if not model:
        raise SystemExit("no membrane_model.json: train the model first (python -m mpstruc.membrane)")
    cands = data.get("candidates") or []
    got = {pdb: (p, why, table) for pdb, p, why, table in _scored([c["pdb"] for c in cands], None, model)}
    n = {"candidates": len(cands), "with_table": 0, "without_table": 0, "why_changed": 0, "why_kept_score_moved": 0}
    for c in cands:
        hit = got.get(c["pdb"])
        if not hit or not hit[2]:
            n["without_table"] += 1
            c.pop("tmEvidence", None)
            continue
        p, why, table = hit
        c["tmEvidence"] = table
        n["with_table"] += 1
        if why != c.get("membraneWhy"):
            if abs(p - (c.get("membrane") or 0)) < 0.0015:
                c["membraneWhy"] = why
                n["why_changed"] += 1
            else:  # annotations moved since the weekly run: the stored score and text stay together
                n["why_kept_score_moved"] += 1
    # The PDB and UniProt records move; the site says when the tables were last checked
    # against them (a weekly run's tables are as old as "generated").
    data["tmChecked"] = (today or datetime.date.today()).isoformat()
    path.write_text(json.dumps(data, separators=(",", ":")))
    return n


def report_md(r: dict, model: dict) -> str:
    t = r["test"]
    L = ["# How likely is a structure to be a membrane protein?", "",
         "A structure counts as a membrane protein if mpstruc lists it, if PDBTM or MemProtMD (which screened every "
         "PDB entry for transmembrane proteins) do, or if OPM classes it as transmembrane. Structures OPM classes only as "
         "monotopic/peripheral or peptide, with no other source, are left out. Signals include UniProt keywords and "
         "families (including the beta-barrel and outer-membrane keywords) and TCDB classes. The strongest signal is "
         "whether the part that crosses the membrane is resolved in the structure: UniProt transmembrane segments, or "
         "for beta-barrels a UniProt barrel domain (such as 'TBDR beta-barrel'). The list of barrel domain names was "
         "chosen from the training period and from published domain definitions, without looking at the test period.", "",
         "Facts from the structure itself (what is resolved, detergents, lipids, the title) and from mpstruc are as they "
         "were on the release day. UniProt, TCDB and InterPro are read as they are today, and some of their annotations "
         "may have been added after (or because of) the structure. So the test-period numbers below may look somewhat "
         "better than they will be for brand-new entries.", "",
         f"It was trained on PDB releases {TRAIN[0]} to {TRAIN[1]} (calibrated on structures each fold had not seen) and "
         f"tested on {TEST[0]} to {TEST[1]}.", "",
         "## Test period", "",
         f"- Structures: {t['n']:,}, of which membrane proteins: {t['membrane']:,}",
         f"- Ranking quality (AUC): {t['auc']}; Brier score: {t['brier']}",
         f"- Called membrane protein (≥50%): {t['at_50']['called_membrane']:,}; of those, {t['at_50']['precision_pct']}% are; "
         f"they include {t['at_50']['recall_pct']}% of all membrane proteins",
         f"- Of the {r['mpstruc_later']['n']} structures mpstruc added from Jul–Dec 2025 releases, "
         f"{r['mpstruc_later']['at_least_50_pct']}% score at least 50% and {r['mpstruc_later']['at_least_10_pct']}% at least 10%",
         ""]
    def kind_name(k):
        return k[0] + k[1:].lower() if k.isupper() else k
    if t.get("by_kind"):
        L += ["### Which kinds are found?", "",
              "Share of each kind in the test period that scores at least 50%: every structure in each mpstruc "
              "group, and the membrane proteins UniProt places in a bacterial outer membrane.", "",
              "| Kind | Structures | Scored 50% or more |", "|---|---|---|"]
        L += [f"| {kind_name(k)} | {v['n']:,} | {v['recall_pct']}% |" for k, v in t["by_kind"].items()]
        L.append("")
    fp = t.get("false_positives")
    if fp:
        L += ["### Wrong calls: model errors or missing labels?", "",
              f"{fp['total']} test structures labelled 'not a membrane protein' score 50% or more. "
              f"{fp['outer_membrane_barrel']} of them are proteins UniProt calls outer-membrane beta-barrels (cell outer "
              "membrane and transmembrane beta strand), which the label sources (mpstruc, PDBTM, MemProtMD, OPM) "
              "probably missed; the same protein is sometimes labelled a membrane protein in another structure. "
              f"The other {fp['other']} are real mistakes. The labels were not changed, so the precision above counts "
              "both kinds as wrong."
              + (f" Outer-membrane barrels counted as wrong: {', '.join(fp['outer_membrane_barrel_ids'])}."
                 if fp["outer_membrane_barrel_ids"] else ""), ""]
    rule = r.get("tm_rule")
    if rule:
        tr, te = rule["train"], rule["test"]
        L += ["### Resolved membrane-crossing segments (the curators' rule)", "",
              "Curators count a structure as a membrane protein only when a significant membrane-crossing part is "
              "resolved in it. How often structures were labelled membrane proteins, by what is resolved: the most "
              f"resolved UniProt transmembrane segments of at least {CORE_MIN_LENGTH} residues in one protein UniProt marks "
              f"transmembrane, in a main protein (at least {round(100 * CORE_MIN_SHARE)}% of the residues modelled) or only "
              "in a minor partner.", "",
              "| Resolved | Training structures | Membrane proteins | Test structures | Membrane proteins |", "|---|---|---|---|---|"]
        for k in sorted(set(tr["segments"]) | set(te["segments"])):
            a, b = tr["segments"].get(k), te["segments"].get(k)
            cell = lambda v: f"{v['n']:,} | {v['membrane_pct']}%" if v else "0 | –"
            L.append(f"| {k} | {cell(a)} | {cell(b)} |")
        L += ["", f"The rule, chosen on the training period only: {CORE_MIN_SEGMENTS} or more such segments in a main "
              f"protein (counted in one chain), or a resolved membrane beta-barrel, make a structure at least "
              f"{round(100 * rule['floor'])}% likely, whatever its protein family or UniProt's 'secreted' or 'cytoplasm' "
              "keywords say. The rule is a floor only, not a weight the model learns, so the model's other weights are "
              "unchanged by it. The floor comes from the training period's held-out predictions (the share of structures "
              "the rule lifts that are labelled membrane proteins); the test period below confirms it but was not used to "
              "choose it. One or two segments are not enough: pili, amyloid fibrils, phage coats, "
              "viral envelope anchors and helices UniProt marks in soluble enzymes have them. "
              f"In the training period (held-out predictions) the rule took {tr['lifted']} structures from below 50% to "
              f"{round(100 * rule['floor'])}%, of which {tr['lifted_membrane_pct']}% are labelled membrane proteins; in the "
              f"test period {te['lifted']}, of which {te['lifted_membrane_pct']}% are. It also raised {tr['raised'] - tr['lifted']} "
              f"training and {te['raised'] - te['lifted']} test structures that were already between 50% and "
              f"{round(100 * rule['floor'])}% (not counted as lifted). Called at 50% without the rule: "
              f"{t['model_alone_at_50']['called_membrane']:,} ({t['model_alone_at_50']['precision_pct']}% precision, "
              f"{t['model_alone_at_50']['recall_pct']}% recall); with it, see above. Lifted in the test period but "
              "labelled 'not a membrane protein': " + (", ".join(te["lifted_not_labelled"]) or "none") + ".", ""]
    misses, spot = t.get("beta_barrel_misses"), r.get("spotlight")
    if spot:
        L += ["### Structures curators asked about", "",
              "| PDB | Shown | Why |", "|---|---|---|"]
        L += [f"| {pdb} | {round(100 * v['p'])}% | {'; '.join(v['why'])} |" for pdb, v in spot.items()]
        L.append("")
    if misses is not None:
        L += ["### Beta-barrels", ""]
        if misses:
            L += [f"mpstruc beta-barrel structures in the test period that still score below 50% ({len(misses)}). "
                  "Barrels whose UniProt entry names no barrel domain (only a plug, or nothing) are still often missed:", "",
                  "| PDB | Shown | Title |", "|---|---|---|"]
            L += [f"| {m['pdb']} | {round(100 * m['p'])}% | {m['title']} |" for m in misses]
            L.append("")
    trade = r.get("monotopic_tradeoff")
    if trade:
        L += ["### Monotopic proteins: membrane proteins or not?", "",
              "mpstruc also lists monotopic proteins, which sit on the membrane without crossing it. If a membrane "
              "protein must have a part that crosses the membrane, they should not count. The same model was trained "
              "three ways; each is tested on its own labels, so the overall numbers differ partly because the "
              f"question differs. This model uses: {r.get('monotopic_label', MONOTOPIC_LABEL)}.", "",
              "| Monotopic structures labelled | AUC | Precision at 50% | Recall at 50% | Monotopic ≥50% | "
              "Alpha-helical ≥50% | Beta-barrel ≥50% |", "|---|---|---|---|---|---|---|"]
        for mode, v in trade.items():
            bk = v["by_kind"]
            g = lambda k: f"{bk[k]['recall_pct']}%" if k in bk else "–"
            L.append(f"| {mode} | {v['auc']} | {v['at_50']['precision_pct']}% | {v['at_50']['recall_pct']}% | "
                     f"{g(MONOTOPIC)} | {g('TRANSMEMBRANE PROTEINS: ALPHA-HELICAL')} | {g('TRANSMEMBRANE PROTEINS: BETA-BARREL')} |")
        L.append("")
    L += ["### Does the percentage mean what it says?", "",
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
    import sys
    if sys.argv[1:2] == ["table"]:  # python -m mpstruc.membrane table: backfill tmEvidence
        print(json.dumps(backfill_table(Path(sys.argv[2]) if len(sys.argv) > 2 else None)))
        raise SystemExit(0)
    # --tradeoff: also train with the other monotopic labellings and report them side by side.
    other = [m for m in ("positive", "exclude", "negative") if m != MONOTOPIC_LABEL] if "--tradeoff" in sys.argv else []
    trade = {m: _brief(train_and_test(monotopic=m, write=False)) for m in other} or None
    res = train_and_test(tradeoff=trade)
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != 'calibration'} for k, v in res.items() if k != "examples" and isinstance(v, dict) and k != "spotlight"}, indent=1))
    print(json.dumps(res["test"]["calibration"], indent=1))
