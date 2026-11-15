"""Count the strands of a membrane beta-barrel in a structure, from its coordinates.

    python -m mpstruc.barrel validate     # expected vs counted on official mpstruc entries -> reports/barrel_validation.md
    python -m mpstruc.barrel 4C00 8US1    # print the count for some structures

Curators place barrel proteins by the number of strands ("Outer Membrane
Autotransporters: typically 12-stranded; Omp85-TpsB: typically 16-stranded",
see curator_rules.json), and UniProt often names no barrel domain for an Omp85
protein such as BamA, so the barrel has to be found in the structure itself.

How a barrel is counted. Only the C-alpha atoms are used (fetched per entity
from the PDB's ModelServer and cached), so it works at any resolution:
 1. Two residues are bridged when their C-alphas and those of their neighbours
    lie side by side (antiparallel: i-1/j+1, i/j, i+1/j-1; parallel: i-1/j-1,
    i/j, i+1/j+1) and both backbones are extended. Runs of bridged residues are
    strands; two runs a beta-bulge splits are merged when they continue each
    other and lean on the same neighbour.
 2. Strands are nodes of a graph, with an edge where two strands are bridged.
    A barrel is a ring in that graph: every strand hydrogen-bonded to the next
    and the last to the first. The count is the number of strands in the
    longest ring made by one chain. Sheets that are not part of the ring
    (POTRA domains, the plug of a TonB-dependent receptor, periplasmic and
    passenger domains) are in other parts of the graph and are not counted.
 3. A barrel whose seam is open (BamA's lateral gate, a barrel caught while it
    folds) has no ring; then the longest row of strands is taken if it wraps
    at least OPEN_MIN_WRAP degrees around an axis at an even radius, and is
    reported as "open". The count of an open barrel is a lower bound ("at
    least N"): strands at the seam may be missing from the row.
 4. Barrels made by several chains (trimeric autotransporters: 4 strands from
    each of 3 chains; TolC, MspA, secretins) are found the same way over the
    copies of one protein, in the deposited model and in the biological
    assembly when that has more copies (the larger ring is kept), and reported
    with the chains they span and where the ring is in the protein.
 5. "upTo": a closed ring that two or more further strands could be part of
    (they bridge two neighbours of the ring more strongly than those are
    bridged to each other) is reported with both numbers; which is the
    barrel is then not decided here.

The count is geometry only. Whether it is a *membrane* barrel is decided by the
caller from UniProt (membrane.counted_barrels: outer membrane, transmembrane
beta strand, Omp85 family, and the ring must lie in the protein's membrane
anchor where UniProt names one): GFP, lipocalins, streptavidin and the head
domains of trimeric autotransporters are closed barrels too.
"""
from __future__ import annotations

import gzip
import json
import math
import re
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

from . import rcsb

VERSION = 6           # bump when the counting changes: cached counts of older versions are recomputed
MODEL_SERVER = "https://models.rcsb.org/v1/{pdb}/atoms?label_atom_id=CA&model_nums=1&encoding=cif"
MAX_BYTES = 30_000_000    # a download larger than this is not a single barrel protein (ribosomes, capsids)
MAX_RESIDUES = 60_000

BRIDGE = 5.9          # C-alpha distance of bridged residues (4.4-5.6 in a regular sheet), in Angstrom
EXTENDED = 5.6        # C-alpha(i-1) to C-alpha(i+1) in an extended backbone (6.4-6.9; a helix has 5.4)
MIN_BRIDGE_PAIRS = 3  # residue pairs between two strands for them to count as neighbours
CHORD_PAIRS = 6       # ... for them to count as well bridged (a ring or row of strands has no such chord)
MIN_STRAND = 4        # residues in a strand (after merging bulges)
MAX_BULGE = 4         # unbridged residues between two runs that are one strand
HAIRPIN_SPAN = 400    # residues apart (in one unbroken run) two bridged residues of a hairpin can be
OPEN_MIN_STRANDS = 8  # an open (unclosed) row of strands shorter than this is not reported
OPEN_MIN_WRAP = 220.0  # degrees it must go around its axis
OPEN_MAX_SPREAD = 0.45  # ... at an even radius: spread (sd / mean) of the strands' distance from the axis
# End strands of an open row nearer to or further from the axis than this x the median are
# dropped. Only strands inside the barrel are (loop 6 of BamA): the first strands of an
# open Omp85 barrel stand away from the wall (BamA or Sam50 with a substrate), and dropping
# them undercounted 9 of 11 such barrels by up to 4 strands (reports/barrel_validation.md).
TRIM_RADIUS = (0.7, math.inf)
DETOUR_MAX = 2        # strands in a detour between two neighbours of a ring (see ring_detours)
SEARCH_STEPS = 400_000  # work limit for the ring search (barrel graphs are nearly rings; this is never reached)

METHOD = "C-alpha geometry"
CARD_FIELDS = ("chain", "count", "method", "closed", "uniprot", "chains", "perChain", "upTo", "from")  # of a barrel, kept on a candidate


def strands_text(b: dict) -> str:
    """'16', 'at least 13' (seam open: strands may be missing from the count) or '8 to 10'
    (a larger ring fits too), for a counted barrel."""
    if not b.get("closed", True):
        return f"at least {b['count']}"
    if b.get("upTo", 0) > b["count"]:
        return f"{b['count']} to {b['upTo']}"
    return str(b["count"])


def barrel_text(b: dict) -> str:
    """'a 16-strand beta-barrel', 'a beta-barrel of at least 13 strands', 'a beta-barrel of 8 to 10 strands'."""
    n = b["count"]
    if not b.get("closed", True) or b.get("upTo", 0) > n:
        return f"a beta-barrel of {strands_text(b)} strands"
    return f"{'an' if str(n)[0] == '8' or n in (11, 18) else 'a'} {n}-strand beta-barrel"


# ---------------------------------------------------------------- coordinates

def _tokens(line: str) -> list[str]:
    if "'" not in line and '"' not in line:
        return line.split()
    return [a or b or c for a, b, c in re.findall(r"'((?:[^']|'(?=\S))*)'|\"((?:[^\"]|\"(?=\S))*)\"|(\S+)", line)]


def parse_cif(text: str, wanted: tuple[str, ...]) -> dict[str, dict[str, list[str]]]:
    """{category: {item: [values]}} for the wanted categories of an mmCIF file (enough of
    the format for what the PDB writes: loops, single values, quoted and ';' text values)."""
    out: dict[str, dict[str, list[str]]] = {}
    lines = text.splitlines()
    i, n = 0, len(lines)
    cols: list[tuple[str, str]] = []
    in_loop = False
    pending: list[str] = []

    def flush_row(vals: list[str]) -> None:
        for (cat, item), v in zip(cols, vals):
            out[cat][item].append(v)

    while i < n:
        line = lines[i]
        if line.startswith(";"):  # multi-line text value
            text_val = [line[1:]]
            i += 1
            while i < n and not lines[i].startswith(";"):
                text_val.append(lines[i])
                i += 1
            i += 1
            pending.append("\n".join(text_val))
        elif line.startswith("#") or not line.strip():
            in_loop, cols, pending = False, [], []
            i += 1
            continue
        elif line.startswith("loop_"):
            in_loop, cols, pending = True, [], []
            i += 1
            continue
        elif line.startswith("_"):
            parts = _tokens(line)
            cat, _, item = parts[0][1:].partition(".")
            keep = cat in wanted
            if in_loop:
                cols.append((cat, item))
                if keep:
                    out.setdefault(cat, {}).setdefault(item, [])
            else:
                cols = [(cat, item)]
                if keep:
                    out.setdefault(cat, {}).setdefault(item, [])
                pending = parts[1:]
                if not pending:
                    i += 1
                    continue
            i += 1
            if in_loop:
                continue
        else:
            if not cols or cols[0][0] not in wanted:
                i += 1
                continue
            pending += _tokens(line)
            i += 1
        if cols and cols[0][0] in wanted:
            while len(pending) >= len(cols):
                flush_row(pending[:len(cols)])
                pending = pending[len(cols):]
        else:
            pending = []
    return out


def _download(url: str) -> bytes | None:
    req = urllib.request.Request(url, headers={"User-Agent": "mpstruc-curation-tool"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                data = r.read(MAX_BYTES + 1)
            return None if len(data) > MAX_BYTES else data
        except urllib.error.HTTPError as err:
            if err.code in (400, 404):
                return None
            if attempt == 4:
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError):
            if attempt == 4:
                raise
        time.sleep(2 ** (attempt + 1))
    return None


def fetch_ca(pdb: str, entities: list[str] | None = None) -> dict | None:
    """C-alpha coordinates of model 1: {"chains": {chain: {"entity": id, "seq": [numbers in
    the entity's sequence], "xyz": [[x, y, z], ...]}}, "operators": [[[3x3 rotation], [shift]],
    ...] of assembly 1 (identity only when the deposited model is the assembly)}.

    entities: only these polymer entities (so one protein of a huge complex is a small
    download). The raw files are cached under cache/barrel_ca/. None if the PDB has no
    such entry or it is too large."""
    pdb = pdb.upper()
    cache = rcsb.CACHE_DIR / "barrel_ca"
    chains: dict[str, dict] = {}
    operators = None
    for ent in (sorted(set(entities)) if entities else [None]):
        path = cache / (f"{pdb}_{ent}.cif.gz" if ent else f"{pdb}.cif.gz")
        if path.exists():
            raw = gzip.decompress(path.read_bytes())
        else:
            raw = _download(MODEL_SERVER.format(pdb=pdb.lower()) + (f"&label_entity_id={ent}" if ent else ""))
            cache.mkdir(parents=True, exist_ok=True)
            path.write_bytes(gzip.compress(raw or b""))
        if not raw:
            if ent is None:
                return None
            continue
        cif = parse_cif(raw.decode("utf-8", "replace"),
                        ("atom_site", "pdbx_struct_assembly_gen", "pdbx_struct_oper_list"))
        a = cif.get("atom_site") or {}
        seen = set()
        for k in range(len(a.get("group_PDB", []))):
            if a["group_PDB"][k] != "ATOM" or a["label_atom_id"][k] != "CA" or not a["label_seq_id"][k].isdigit():
                continue
            ch, seq = a["auth_asym_id"][k], int(a["label_seq_id"][k])
            if (ch, seq) in seen:  # alternate conformations: the first one
                continue
            seen.add((ch, seq))
            c = chains.setdefault(ch, {"entity": a["label_entity_id"][k], "label": a["label_asym_id"][k], "seq": [], "xyz": []})
            c["seq"].append(seq)
            c["xyz"].append([float(a["Cartn_x"][k]), float(a["Cartn_y"][k]), float(a["Cartn_z"][k])])
        if operators is None:
            operators = _assembly_operators(cif, {c["label"] for c in chains.values()})
    if not chains:
        return None
    return {"chains": chains, "operators": operators or []}


def _assembly_operators(cif: dict, labels: set[str]) -> list:
    """The operators of the first assembly that are applied to these chains, when it is
    built from more than the identity (a trimer on a crystal axis, say)."""
    gen, ops = cif.get("pdbx_struct_assembly_gen") or {}, cif.get("pdbx_struct_oper_list") or {}
    if not gen or not ops:
        return []
    by_id = {}
    for k, oid in enumerate(ops.get("id", [])):
        try:
            rot = [[float(ops[f"matrix[{r}][{c}]"][k]) for c in (1, 2, 3)] for r in (1, 2, 3)]
            vec = [float(ops[f"vector[{r}]"][k]) for r in (1, 2, 3)]
        except (KeyError, ValueError):
            continue
        by_id[oid] = [rot, vec]
    first = (gen.get("assembly_id") or ["1"])[0]
    out = []
    for k, aid in enumerate(gen.get("assembly_id", [])):
        if aid != first or not labels & set(gen["asym_id_list"][k].split(",")):
            continue
        expr = gen["oper_expression"][k]
        if ")(" in expr:  # products of operators (capsids): not a barrel protein's assembly
            return []
        for part in expr.strip("()").split(","):
            lo, _, hi = part.partition("-")
            ids = [str(x) for x in range(int(lo), int(hi) + 1)] if hi and lo.isdigit() and hi.isdigit() else [part]
            out += [by_id[x] for x in ids if x in by_id and by_id[x] not in out]
    identity = [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [0.0, 0.0, 0.0]]
    return [] if out in ([], [identity]) else out


# ---------------------------------------------------------------- strands

def _near_pairs(xyz, cut: float):
    """(i, j) with i < j and |xyz[i] - xyz[j]| < cut, by a cell grid."""
    import numpy as np
    n = len(xyz)
    cell = np.floor(xyz / cut).astype(np.int64)
    cell -= cell.min(0) - 1
    dims = cell.max(0) + 2
    key = (cell[:, 0] * dims[1] + cell[:, 1]) * dims[2] + cell[:, 2]
    order = np.argsort(key, kind="stable")
    skey = key[order]
    ii, jj = [], []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                nk = key + (dx * dims[1] + dy) * dims[2] + dz
                lo, hi = np.searchsorted(skey, nk, "left"), np.searchsorted(skey, nk, "right")
                cnt = hi - lo
                if not cnt.sum():
                    continue
                a = np.repeat(np.arange(n), cnt)
                within = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
                b = order[np.repeat(lo, cnt) + within]
                m = a < b
                a, b = a[m], b[m]
                d = np.linalg.norm(xyz[a] - xyz[b], axis=1)
                m = d < cut
                ii.append(a[m])
                jj.append(b[m])
    if not ii:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    return np.concatenate(ii), np.concatenate(jj)


def find_strands(chains: dict[str, dict]) -> tuple[list[dict], dict[tuple[int, int], int]]:
    """Strands and which of them are bridged.

    chains: {chain: {"seq": [residue numbers], "xyz": [[x, y, z], ...]}} (C-alpha atoms).
    Returns ([{"chain", "first", "last", "n" (bridged residues)}, ...] in chain and
    sequence order, {(a, b): residue pairs bridged between strands a < b})."""
    import numpy as np
    names = list(chains)
    if not names:
        return [], {}
    xyz = np.array([p for c in names for p in chains[c]["xyz"]], dtype=float).reshape(-1, 3)
    seq = np.array([s for c in names for s in chains[c]["seq"]], dtype=np.int64)
    cid = np.array([k for k, c in enumerate(names) for _ in chains[c]["seq"]], dtype=np.int64)
    n = len(seq)
    if n < 6:
        return [], {}
    idx = np.arange(n)
    has_next = np.zeros(n, dtype=bool)
    has_next[:-1] = (cid[1:] == cid[:-1]) & (seq[1:] == seq[:-1] + 1)
    has_prev = np.zeros(n, dtype=bool)
    has_prev[1:] = has_next[:-1]
    inner = has_next & has_prev
    ext = np.zeros(n, dtype=bool)
    mid = idx[inner]
    ext[mid] = np.linalg.norm(xyz[mid + 1] - xyz[mid - 1], axis=1) > EXTENDED
    I, J = _near_pairs(xyz, BRIDGE)
    keys = np.sort(I * n + J)

    def near(a, b):
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        k = lo * n + hi
        pos = np.searchsorted(keys, k)
        pos[pos >= len(keys)] = len(keys) - 1
        return (keys[pos] == k) if len(keys) else np.zeros(len(k), dtype=bool)

    ok = ext[I] & ext[J]
    I, J = I[ok], J[ok]
    same = cid[I] == cid[J]
    # Antiparallel: (i-1, j+1), (i, j), (i+1, j-1); in one chain the two strands need a turn between them.
    anti = near(I - 1, J + 1) & near(I + 1, J - 1) & (~same | (J - I >= 4))
    par = near(I - 1, J - 1) & near(I + 1, J + 1) & (~same | (J - I >= 4))
    pairs = set()
    for a, b in zip(I[anti].tolist(), J[anti].tolist()):
        pairs.update(((a - 1, b + 1), (a, b), (a + 1, b - 1)))
    for a, b in zip(I[par].tolist(), J[par].tolist()):
        pairs.update(((a - 1, b - 1), (a, b), (a + 1, b + 1)))
    marked = sorted({r for p in pairs for r in p})
    # Runs of bridged residues.
    raw_runs: list[list[int]] = []
    for r in marked:
        if raw_runs and raw_runs[-1][-1] == r - 1 and has_prev[r]:
            raw_runs[-1].append(r)
        else:
            raw_runs.append([r])
    # A run whose residues are bridged to each other is a hairpin with a tight turn (its
    # turn residues touch a neighbouring turn): two strands, split at the turn.
    start_of = {r: k for k, run in enumerate(raw_runs) for r in run}
    turns: dict[int, list[float]] = defaultdict(list)
    for a, b in pairs:
        if start_of[a] == start_of[b] and abs(a - b) <= HAIRPIN_SPAN:
            turns[start_of[a]].append((a + b) / 2)
    run_of: dict[int, int] = {}
    runs: list[list[int]] = []
    for k, run in enumerate(raw_runs):
        cuts: list[float] = []
        for t in sorted(turns.get(k, ())):  # pairs of one hairpin share its turn's centre
            if not cuts or t - cuts[-1] > 2:
                cuts.append(t)
        piece: list[int] = []
        for r in run:
            if cuts and r > cuts[0]:
                cuts.pop(0)
                if piece:
                    runs.append(piece)
                piece = []
            piece.append(r)
        if piece:
            runs.append(piece)
    for k, run in enumerate(runs):
        for r in run:
            run_of[r] = k
    link: dict[tuple[int, int], set] = defaultdict(set)
    for a, b in pairs:
        ra, rb = run_of[a], run_of[b]
        if ra != rb:
            link[(min(ra, rb), max(ra, rb))].add((a, b))
    partners: dict[int, set[int]] = defaultdict(set)
    for (a, b), v in link.items():
        if len(v) >= 2:
            partners[a].add(b)
            partners[b].add(a)
    # Merge runs a beta-bulge split: next to each other in the chain, not bridged to each
    # other (that is a hairpin), and bridged to a common neighbour.
    parent = list(range(len(runs)))

    def root(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for k in range(len(runs) - 1):
        a, b = runs[k], runs[k + 1]
        if cid[a[-1]] != cid[b[0]] or seq[b[0]] - seq[a[-1]] - 1 > MAX_BULGE:
            continue
        if (k, k + 1) in link or not partners[k] & partners[k + 1]:
            continue
        parent[root(k + 1)] = root(k)
    groups: dict[int, list[int]] = defaultdict(list)
    for k in range(len(runs)):
        groups[root(k)].append(k)
    strands, new_id = [], {}
    for g in sorted(groups):
        res = [r for k in groups[g] for r in runs[k]]
        if len(res) < MIN_STRAND:
            continue
        for k in groups[g]:
            new_id[k] = len(strands)
        strands.append({"chain": names[cid[res[0]]], "first": int(seq[res[0]]), "last": int(seq[res[-1]]), "n": len(res),
                        "centre": xyz[res].mean(0).tolist()})
    edges: dict[tuple[int, int], int] = defaultdict(int)
    for (a, b), v in link.items():
        if a in new_id and b in new_id and new_id[a] != new_id[b]:
            x, y = sorted((new_id[a], new_id[b]))
            edges[(x, y)] += len(v)
    return strands, {e: k for e, k in edges.items() if k >= MIN_BRIDGE_PAIRS}


# ---------------------------------------------------------------- rings

def _adjacency(nodes: set[int], edges: dict[tuple[int, int], int]) -> dict[int, set[int]]:
    adj: dict[int, set[int]] = {v: set() for v in nodes}
    for a, b in edges:
        if a in nodes and b in nodes:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def _strong(nodes: set[int], edges: dict[tuple[int, int], int]) -> dict[int, set[int]]:
    return _adjacency(nodes, {e: k for e, k in edges.items() if k >= CHORD_PAIRS})


def longest_ring(nodes: set[int], edges: dict[tuple[int, int], int], meander: bool = True) -> list[int]:
    """The longest ring of at least 4 strands in the strand graph, [] if none.

    In a barrel a strand is bridged to the one before and the one after only, so the ring
    must have no chord: no well-bridged pair (CHORD_PAIRS residue pairs or more) other than
    neighbours in the ring. A short strand lying across two neighbours (the C-terminal tail
    of invasin) or a hairpin leaning on the wall is therefore left out.

    meander: the strands follow each other around the ring in the order they have in the
    chain (strand numbers are in chain order), as in every outer-membrane barrel; a hairpin
    folded back into the lumen (loop 6 of BamA) is then never part of it. Off for barrels
    made by several chains."""
    adj = _adjacency(nodes, edges)
    while True:  # strands with fewer than two neighbours cannot be in a ring
        drop = [v for v, nb in adj.items() if len(nb) < 2]
        if not drop:
            break
        for v in drop:
            for u in adj.pop(v):
                if u in adj:
                    adj[u].discard(v)
    strong = _strong(set(adj), edges)
    best: list[int] = []
    steps = 0
    for start in sorted(adj):
        path, on = [start], {start}
        stack = [iter(sorted(u for u in adj[start] if u > start))]
        while stack:
            steps += 1
            if steps > SEARCH_STEPS:
                return best
            u = next(stack[-1], None)
            if u is None:
                stack.pop()
                on.discard(path.pop())
                continue
            if u in on or (strong[u] & on) - {path[-1], start} or (meander and u < path[-1]):
                continue
            closes = start in adj[u] and len(path) >= 3
            if start in strong[u] and len(path) > 1 and not closes:
                continue  # touches the first strand too early: a chord
            path.append(u)
            on.add(u)
            if closes and len(path) > len(best):
                best = list(path)
            if start in strong[u] and len(path) > 2:  # the ring is closed here; going on would make this a chord
                on.discard(path.pop())
                continue
            stack.append(iter(sorted(w for w in adj[u] if w > start)))
    return best


def longest_row(nodes: set[int], edges: dict[tuple[int, int], int], meander: bool = True) -> list[int]:
    """The longest row of strands, each bridged to the next and to no other strand of the
    row (see longest_ring)."""
    adj = _adjacency(nodes, edges)
    strong = _strong(nodes, edges)
    best: list[int] = []
    steps = 0
    for start in sorted(adj):
        path, on = [start], {start}
        stack = [iter(sorted(adj[start]))]
        while stack:
            steps += 1
            if steps > SEARCH_STEPS:
                return best
            u = next(stack[-1], None)
            if u is None:
                if len(path) > len(best):
                    best = list(path)
                stack.pop()
                on.discard(path.pop())
                continue
            if u in on or (strong[u] & on) - {path[-1]} or (meander and u < path[-1]):
                continue
            path.append(u)
            on.add(u)
            stack.append(iter(sorted(adj[u])))
    return best


def wrap(centres: list[list[float]]) -> tuple[float, float]:
    """(degrees a row of strands goes around its axis, spread of their distance from it).

    The axis is the normal of the plane through the strand centres (they lie on a ring in
    a barrel). A flat or twisted sheet gives well under 220 degrees; a sandwich of two
    sheets gives an uneven radius."""
    import numpy as np
    pts = np.array(centres, dtype=float)
    if len(pts) < 4:
        return 0.0, 1.0
    c = pts - pts.mean(0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    u, v = vt[0], vt[1]
    x, y = c @ u, c @ v
    ang = np.arctan2(y, x)
    d = np.diff(ang)
    d = (d + math.pi) % (2 * math.pi) - math.pi
    r = np.hypot(x, y)
    return abs(float(np.degrees(d.sum()))), float(r.std() / r.mean()) if r.mean() else 1.0


def trim_row(row: list[int], centres: list[list[float]]) -> list[int]:
    """A row of strands without end strands that are off the barrel wall: a hairpin folded
    into the lumen (loop 6 of BamA) or a domain outside it leans on the last strand of an
    open barrel, much nearer to or further from the axis than the wall is."""
    import numpy as np
    row, pts = list(row), [list(c) for c in centres]
    while len(row) > 4:
        p = np.array(pts, dtype=float)
        c = p - p.mean(0)
        _, _, vt = np.linalg.svd(c, full_matrices=False)
        r = np.hypot(c @ vt[0], c @ vt[1])
        mid = float(np.median(r))
        off = [k for k in (0, len(row) - 1) if not TRIM_RADIUS[0] * mid <= r[k] <= TRIM_RADIUS[1] * mid]
        if not off:
            break
        k = max(off, key=lambda k: abs(r[k] - mid))
        row.pop(k)
        pts.pop(k)
    return row


def ring_detours(ring: list[int], nodes: set[int], edges: dict[tuple[int, int], int]) -> int:
    """How many more strands a closed ring could have: strands outside it (at most
    DETOUR_MAX in a row, in chain order) that join two neighbours of the ring by bridges
    all stronger than the bridge between those two. Then two rings fit the same strands
    (Chlamydia MOMP 9OQ7: 8, or 10 with a hairpin between the last and the first strand)
    and the count is reported as a range."""
    w = lambda a, b: edges.get((min(a, b), max(a, b)), 0)
    on, extra = set(ring), 0
    lo, hi = min(ring), max(ring)
    for i, u in enumerate(ring):
        v = ring[(i + 1) % len(ring)]
        base = w(u, v)
        closing = i == len(ring) - 1
        ok = (lambda x: x > hi or x < lo) if closing else (lambda x: min(u, v) < x < max(u, v))
        free = [x for x in nodes if x not in on and ok(x)]
        found = 0
        for x in free:
            if w(u, x) > base and w(x, v) > base:
                found = 1
                break
        if not found and DETOUR_MAX >= 2:
            for x in free:
                if w(u, x) <= base:
                    continue
                if any(y != x and w(x, y) > base and w(y, v) > base for y in free):
                    found = 2
                    break
        extra += found
    return extra


def _apply(op: list, xyz: list) -> list:
    rot, vec = op
    return [[sum(rot[r][k] * p[k] for k in range(3)) + vec[r] for r in range(3)] for p in xyz]


def count_barrels(ca: dict) -> dict:
    """Barrels in a structure, from fetch_ca output.

    {"chains": {chain: {"entity", "count", "closed", "strands": [[first, last], ...] in the
                entity's numbering, around the barrel[, "upTo": see ring_detours]}} for each
                chain that makes a barrel on its own (a ring; or an open row that wraps
                around, "closed": False),
     "oligomers": [{"entity", "chains": [...], "count", "perChain", "first", "last" (the
                ring's residues in the entity's numbering)}] barrels made by several copies
                of one protein,
     "residues": C-alpha atoms read}."""
    chains = ca["chains"]
    out: dict = {"chains": {}, "oligomers": [], "residues": sum(len(c["seq"]) for c in chains.values())}
    if out["residues"] > MAX_RESIDUES:
        out["skipped"] = "too large"
        return out
    strands, edges = find_strands(chains)
    by_chain: dict[str, set[int]] = defaultdict(set)
    for k, s in enumerate(strands):
        by_chain[s["chain"]].add(k)
    for ch, nodes in by_chain.items():
        ring = longest_ring(nodes, edges)
        closed = True
        if len(ring) < OPEN_MIN_STRANDS:
            # No ring, or only a small one (two hairpins side by side): an open barrel?
            row = longest_row(nodes, edges)
            row = trim_row(row, [strands[k]["centre"] for k in row])
            if len(row) >= max(OPEN_MIN_STRANDS, len(ring) + 1):
                deg, spread = wrap([strands[k]["centre"] for k in row])
                if deg >= OPEN_MIN_WRAP and spread <= OPEN_MAX_SPREAD:
                    ring, closed = row, False
        if ring:
            out["chains"][ch] = {"entity": chains[ch]["entity"], "count": len(ring), "closed": closed,
                                 "strands": [[strands[k]["first"], strands[k]["last"]] for k in ring]}
            more = ring_detours(ring, nodes, edges) if closed else 0
            if more:
                out["chains"][ch]["upTo"] = len(ring) + more
    # Barrels made by several copies of a protein: only looked for in proteins whose
    # chains make no barrel on their own.
    by_entity: dict[str, list[str]] = defaultdict(list)
    for ch, c in chains.items():
        by_entity[c["entity"]].append(ch)
    for ent, members in by_entity.items():
        if any(ch in out["chains"] for ch in members):
            continue
        sub_strands = strands
        nodes = {k for k, s in enumerate(strands) if s["chain"] in members}
        ring = longest_ring(nodes, edges, meander=False) if len(members) > 1 else []
        # The biological assembly has more copies than were deposited (a trimer on a crystal
        # axis; MspA 1UUN, two of its eight chains): the ring is looked for there too, and
        # the larger one kept.
        if ca.get("operators") and len(ca["operators"]) > 1 and \
                sum(len(chains[ch]["seq"]) for ch in members) * len(ca["operators"]) <= MAX_RESIDUES:
            copies = {}
            for n, op in enumerate(ca["operators"]):
                for ch in members:
                    copies[f"{ch}-{n + 1}"] = {**chains[ch], "xyz": _apply(op, chains[ch]["xyz"])}
            all_strands, all_edges = find_strands(copies)
            whole = longest_ring(set(range(len(all_strands))), all_edges, meander=False)
            if len(whole) > len(ring):
                sub_strands, ring = all_strands, whole
        used = sorted({sub_strands[k]["chain"] for k in ring})
        if len(ring) >= 4 and len(used) > 1:
            out["oligomers"].append({"entity": ent, "chains": used, "count": len(ring),
                                     "perChain": round(len(ring) / len(used), 1),
                                     "first": min(sub_strands[k]["first"] for k in ring),
                                     "last": max(sub_strands[k]["last"] for k in ring)})
    return out


# ---------------------------------------------------------------- cached counts

def _cache_path() -> Path:
    return rcsb.CACHE_DIR / "barrel_strands.json"


def barrel_counts(wanted: dict[str, list[str] | None]) -> dict[str, dict | None]:
    """count_barrels for many structures, cached in cache/barrel_strands.json.

    wanted: {PDB ID: polymer entity ids to look at (None: every chain)}. A structure the
    PDB could not serve (or that is too large) maps to None; network errors are not cached."""
    path = _cache_path()
    cache = json.loads(path.read_text()) if path.exists() else {}
    out: dict[str, dict | None] = {}
    dirty = 0
    for pdb, ents in wanted.items():
        key = pdb.upper() + ("" if ents is None else ":" + ",".join(sorted(set(ents))))
        hit = cache.get(key)
        if hit is not None and hit.get("v") == VERSION:
            out[pdb] = hit.get("r")
            continue
        try:
            ca = fetch_ca(pdb, ents)
        except Exception as err:  # the count is simply not known this time
            print(f"barrel: could not fetch {pdb} ({err})", file=sys.stderr)
            out[pdb] = None
            continue
        res = count_barrels(ca) if ca else None
        cache[key] = {"v": VERSION, "r": res}
        out[pdb] = res
        dirty += 1
        if dirty % 25 == 0:
            path.write_text(json.dumps(cache))
    if dirty:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cache))
    return out


def best_barrel(res: dict | None, entities: set[str] | None = None) -> dict | None:
    """The barrel to report for a structure: {"chain", "count", "method", "closed",
    "entity", "first", "last" (its residues in the entity's numbering)[, "upTo"], and for
    a barrel made by several chains also "chains" and "perChain"}.

    entities: only barrels of these polymer entities. A closed single-chain barrel is
    preferred, then the one with the most strands."""
    if not res:
        return None
    cands = []
    for ch, b in sorted((res.get("chains") or {}).items()):
        if entities is None or b["entity"] in entities:
            cands.append({"chain": ch, "count": b["count"], "method": METHOD, "closed": b["closed"], "entity": b["entity"],
                          "first": min(x[0] for x in b["strands"]), "last": max(x[1] for x in b["strands"]),
                          **({"upTo": b["upTo"]} if b.get("upTo") else {})})
    for o in res.get("oligomers") or []:
        if entities is None or o["entity"] in entities:
            cands.append({"chain": o["chains"][0].split("-")[0], "count": o["count"], "method": METHOD, "closed": True,
                          "entity": o["entity"], "chains": len(o["chains"]), "perChain": o["perChain"],
                          "first": o.get("first"), "last": o.get("last")})
    if not cands:
        return None
    return max(cands, key=lambda b: (b["count"] >= 8, b["closed"], b["count"], "chains" not in b))


# ---------------------------------------------------------------- validation

# Strand counts from the papers and from mpstruc's own text, for official mpstruc entries.
# (PDB, protein, expected strands, mpstruc subgroup it is in.)
EXPECTED = [
    # Outer Membrane Autotransporters (G02-S006): 12, trimeric ones 4 per chain x 3
    ("1UYN", "NalP translocator domain", 12, "G02-S006"), ("2QOM", "EspP", 12, "G02-S006"),
    ("3SLJ", "EspP, pre-cleavage", 12, "G02-S006"), ("3KVN", "EstA", 12, "G02-S006"),
    ("3AEH", "Hbp", 12, "G02-S006"), ("3QQ2", "BrkA", 12, "G02-S006"), ("4MEE", "AIDA-I", 12, "G02-S006"),
    ("4E1S", "Intimin", 12, "G02-S006"), ("4E1T", "Invasin", 12, "G02-S006"),
    ("2GR8", "Hia (trimeric, 4 per chain)", 12, "G02-S006"), ("2LME", "YadA (trimeric, 4 per chain)", 12, "G02-S006"),
    ("4C00", "TamA", 16, "G02-S006"), ("6WIL", "CdiB (a TpsB transporter)", 16, "G02-S006"),
    ("6WIM", "CdiB (a TpsB transporter)", 16, "G02-S006"),
    # Omp85-TpsB (G02-S007): 16
    ("4QKY", "FhaC", 16, "G02-S007"), ("3NJT", "FhaC R450A", 16, "G02-S007"), ("4K3B", "BamA", 16, "G02-S007"),
    ("4K3C", "BamA", 16, "G02-S007"), ("4N75", "BamA barrel", 16, "G02-S007"), ("4C4V", "BamA", 16, "G02-S007"),
    ("6QGW", "BamA with nanobody", 16, "G02-S007"), ("5EKQ", "BamACDE", 16, "G02-S007"),
    ("5D0O", "BamABCDE", 16, "G02-S007"), ("5LJO", "BamABCDE (related code)", 16, "G02-S007"),
    ("6V05", "BamABCDE with substrate", 16, "G02-S007"), ("8BWC", "BamABCDE", 16, "G02-S007"),
    ("8US1", "TamA", 16, "G02-S007"),
    # Controls from the other beta-barrel subfamilies
    ("1BXW", "OmpA", 8, "control"), ("1QJP", "OmpA", 8, "control"), ("1THQ", "PagP", 8, "control"),
    ("1QJ8", "OmpX", 8, "control"), ("1K24", "OpcA", 10, "control"), ("1I78", "OmpT", 10, "control"),
    ("1QD6", "OmpLA", 12, "control"), ("1UYN", "NalP", 12, "control"), ("1TLY", "Tsx", 12, "control"),
    ("3DWO", "FadL", 14, "control"), ("1T16", "FadL", 14, "control"), ("2OMF", "OmpF", 16, "control"),
    ("1OPF", "OmpF", 16, "control"), ("2POR", "porin", 16, "control"), ("1PRN", "porin", 16, "control"),
    ("1A0S", "ScrY sucrose porin", 18, "control"), ("1MAL", "LamB maltoporin", 18, "control"),
    ("2O4V", "OprP", 16, "control"), ("3SY7", "OccD1 (OprD)", 18, "control"),
    ("1FEP", "FepA (TonB-dependent)", 22, "control"), ("2FCP", "FhuA (TonB-dependent)", 22, "control"),
    ("1KMO", "FecA (TonB-dependent)", 22, "control"), ("1NQE", "BtuB (TonB-dependent)", 22, "control"),
    ("4Q35", "LptD", 26, "control"), ("5IV9", "LptD", 26, "control"), ("2J1N", "OmpC", 16, "control"),
    ("3EMN", "VDAC1 (mitochondrial)", 19, "control"), ("2JK4", "VDAC1 (mitochondrial)", 19, "control"),
    ("1EK9", "TolC (trimer, 4 per chain)", 12, "control"), ("7AHL", "alpha-hemolysin (heptamer, 2 per chain)", 14, "control"),
]
# Barrels made by several chains where fewer copies were deposited than the assembly has.
EXPECTED += [("1UUN", "MspA (octamer, 2 per chain; two chains deposited)", 16, "control")]
# Omp85 barrels (16 strands) caught with their seam open: a substrate in BamA or Sam50,
# TamA with TamB. Their count is reported as a lower bound; listed apart from the table above.
OPEN_SEAM = [
    ("7RI4", "BamA with a substrate (EspP)", 16), ("7TT3", "BamA with a substrate", 16), ("7TSZ", "BamA with a substrate", 16),
    ("7RJ5", "BamA with a substrate", 16), ("9SOW", "BamA with a substrate (new structure)", 16),
    ("6WUH", "Sam50 with a substrate (Tom40)", 16), ("6WUJ", "Sam50 with a substrate", 16), ("6WUT", "Sam50 with a substrate", 16),
    ("6WUM", "Sam50 with a substrate", 16), ("25WF", "SAM50 (new structure)", 16), ("9TGQ", "TamA with TamB (new structure)", 16),
]
# Closed barrels that are not membrane proteins: counted, but must never be read as one.
SOLUBLE = [("1EMA", "GFP"), ("1STP", "streptavidin"), ("1BBP", "bilin-binding protein (lipocalin)"),
           ("8TIM", "triosephosphate isomerase (TIM barrel)"), ("3PRW", "BamB (beta-propeller)"),
           ("3VE2", "TbpB (lipoprotein on the outer membrane)"), ("3KVD", "fHbp (lipoprotein on the outer membrane)")]
# Soluble fragments of outer-membrane proteins. UniProt's keywords describe the full-length
# protein, so only where the ring is tells them from the membrane barrel
# (membrane.membrane_barrel). None may be read as a membrane barrel.
GUARD = [
    ("1S7M", "Hia binding domain HiaBD1 (trimeric autotransporter head)"), ("3EMI", "Hia 307-422 (trimeric autotransporter)"),
    ("3LA9", "BpaA head domain (trimeric autotransporter)"), ("4USX", "BPSL2063 YadA-like head (trimeric autotransporter)"),
    ("3D9X", "BadA head (trimeric autotransporter)"), ("3WPR", "AtaA stalk fragment fused to GCN4 (trimeric autotransporter)"),
    ("2YO2", "SadA 255-358 fused to GCN4 (trimeric autotransporter)"), ("2YNZ", "SadA 823-947 fused to GCN4 (trimeric autotransporter)"),
    ("2XQH", "EibD immunoglobulin-binding fragment (trimeric autotransporter)"),
    ("2QCZ", "BamA POTRA domains"), ("2QDF", "BamA POTRA domains"), ("3EFC", "BamA periplasmic domain"),
    ("3OG5", "BamA POTRA 4-5"), ("4XGA", "BamB with BamA POTRA 3-5"), ("3MC8", "Omp85 POTRA 1-3 (Anabaena)"),
    ("2X8X", "Omp85 N-terminal domain (Thermosynechococcus)"), ("5UAY", "Toc75 POTRA domains"),
    ("1WXR", "Hbp passenger domain (autotransporter)"), ("1DAB", "pertactin P.69 passenger domain (autotransporter)"),
    ("1F02", "intimin extracellular fragment"),
    ("3BWU", "FimD usher N-terminal domain with FimC-FimH"), ("1ZE3", "FimD usher N-terminal domain with FimC-FimH pilin domain"),
    ("3L48", "PapC usher C-terminal domain"),
    ("1ZZV", "FecA periplasmic signalling domain (TonB-dependent receptor)"),
    ("2A02", "PupA periplasmic signalling domain (TonB-dependent receptor)"),
]


def _as_used(ids: list[str], every=()) -> dict[str, dict]:
    """What the membrane score and the placement rules read in these structures
    (membrane.counted_barrels with every outer-membrane protein counted): {pdb: {accession:
    barrel}}, a structure absent if it could not be counted."""
    from . import membrane
    ext, fam, ents = rcsb.entry_extras(ids), rcsb.entry_families(ids), rcsb.entry_entities(ids)
    return membrane.counted_barrels(ids, ext, fam, None, ents, every=every)


def validate() -> dict:
    ids = sorted({p for p, *_ in EXPECTED} | {p for p, *_ in OPEN_SEAM} | {p for p, _ in SOLUBLE} | {p for p, _ in GUARD})
    got = barrel_counts({p: None for p in ids})
    ruled_ids = {p for p, _, _, sg in EXPECTED if sg != "control"}
    used = _as_used(ids)                       # as the membrane score reads them
    voted = _as_used(sorted(ruled_ids), ruled_ids)  # as the rules read a structure suggested under a ruled subfamily

    def read(pdb: str) -> str:
        """How production reads the structure: the filter the geometry-only count leaves out."""
        if pdb not in used:
            return "not counted"
        if any(b["count"] >= 8 for b in used[pdb].values()):  # membrane.CORE_MIN_STRANDS
            return "yes"
        if any(b.get("vote") for b in (voted.get(pdb) or {}).values()):
            return "only when suggested under a ruled subfamily (UniProt names no outer-membrane protein)"
        return "no"

    def row(pdb, name, want, sg):
        b = best_barrel(got.get(pdb))
        return {"pdb": pdb, "protein": name, "subgroup": sg, "expected": want, "counted": b["count"] if b else None,
                "closed": b["closed"] if b else None, "chains": (b or {}).get("chains", 1),
                "ok": bool(b) and b["count"] == want, "membrane": read(pdb)}

    rows = [row(*x) for x in EXPECTED]
    seen, uniq = set(), []
    for r in rows:  # 1UYN is listed both as an autotransporter and as a control
        if r["pdb"] not in seen:
            seen.add(r["pdb"])
            uniq.append(r)
    open_rows = [row(p, name, want, "open seam") for p, name, want in OPEN_SEAM]
    all_open = [r for r in uniq + open_rows if r["closed"] is False]

    def soluble(pairs):
        out = []
        for pdb, name in pairs:
            b = best_barrel(got.get(pdb))
            out.append({"pdb": pdb, "protein": name, "counted": b["count"] if b else None, "closed": b["closed"] if b else None,
                        "chains": (b or {}).get("chains", 1), "membrane": read(pdb)})
        return out

    sol, guard = soluble(SOLUBLE), soluble(GUARD)
    ok = sum(r["ok"] for r in uniq)
    ruled = [r for r in uniq if r["subgroup"] != "control"]
    return {"rows": uniq, "open": open_rows, "soluble": sol, "guard": guard,
            "n": len(uniq), "correct": ok, "accuracy": round(ok / len(uniq), 3),
            "ruled_n": len(ruled), "ruled_correct": sum(r["ok"] for r in ruled),
            "within_one": sum(1 for r in uniq if r["counted"] is not None and abs(r["counted"] - r["expected"]) <= 1),
            "open_n": len(all_open), "open_correct": sum(r["ok"] for r in all_open),
            "guard_failed": [r["pdb"] for r in sol + guard if r["membrane"] != "no"]}


def validation_md(v: dict) -> str:
    def line(r):
        kind = "-" if r["counted"] is None else ("closed" if r["closed"] else "open seam") + \
            (f", {r['chains']} chains" if r["chains"] > 1 else "")
        return (f"| {r['pdb']} | {r['protein']} | {r['subgroup']} | {r['expected']} | "
                f"{'none found' if r['counted'] is None else r['counted']} | {kind} | {'ok' if r['ok'] else '**differs**'} | {r['membrane']} |")
    head = ["| PDB | Protein | In mpstruc under | Expected | Counted | Barrel | | Read as a membrane barrel |", "|---|---|---|---|---|---|---|---|"]
    L = ["# Counting barrel strands: expected vs counted", "",
         "Strands of the membrane beta-barrel counted from C-alpha coordinates (`mpstruc/barrel.py`), against the",
         "number in the literature, for official mpstruc entries. Regenerate with `python -m mpstruc.barrel validate`.", "",
         f"**Exactly right: {v['correct']} of {v['n']} ({round(100 * v['accuracy'])}%)**; within one strand: "
         f"{v['within_one']} of {v['n']}. In the two subfamilies with a curator rule: {v['ruled_correct']} of {v['ruled_n']}.", "",
         "\"Counted\" is geometry only, over every chain. \"Read as a membrane barrel\" is what the membrane score and the",
         "placement rules do with it (`membrane.counted_barrels`): only chains of proteins UniProt places in the outer",
         "membrane are looked at, and the ring must lie in the protein's membrane part.", ""] + head
    L += [line(r) for r in v["rows"]]
    L += ["", "## Barrels with an open seam", "",
          f"A barrel caught with a substrate in it or while it folds has no closed ring; its count is a lower bound and is",
          f"shown as \"at least N\". Exactly right: **{v['open_correct']} of {v['open_n']}** open barrels (those in the table above",
          f"and the Omp85 barrels below, 16 strands each).", ""] + head
    L += [line(r) for r in v["open"]]

    def table(rows):
        out = ["| PDB | Protein | Counted (geometry) | Read as a membrane barrel |", "|---|---|---|---|"]
        for r in rows:
            what = "no barrel" if r["counted"] is None else \
                f"{r['counted']} ({'closed' if r['closed'] else 'open'}{', ' + str(r['chains']) + ' chains' if r['chains'] > 1 else ''})"
            out.append(f"| {r['pdb']} | {r['protein']} | {what} | {'no' if r['membrane'] == 'no' else '**' + r['membrane'] + '**'} |")
        return out

    L += ["", "## Barrels and sheets that are not membrane barrels", "",
          "The count is geometry only: soluble barrels close too. They are never read as membrane barrels because",
          "the membrane score and the placement rules also require UniProt to place the protein in the outer membrane",
          "(see `membrane.OM_FAMILIES`, `membrane.counted_barrels`).", ""] + table(v["soluble"])
    L += ["", "## Guard set: soluble fragments of outer-membrane proteins", "",
          "UniProt's keywords describe the full-length protein, so a head, stalk, passenger, POTRA or periplasmic",
          "fragment carries them too. A ring of strands in such a fragment is not read as a membrane barrel: it must lie",
          "in the membrane part UniProt names (transmembrane strands, a barrel domain, the membrane anchor or",
          "translocator domain), and a ring made by several chains that cannot be checked that way (UniProt names no",
          "membrane part) needs 12 or more strands, a whole number per chain (`membrane.membrane_barrel`). A ring of fewer",
          "than 8 strands is never a membrane barrel.", "",
          f"**Read as a membrane barrel: {len(v['guard_failed'])} of {len(v['soluble']) + len(v['guard'])}** (must be 0"
          + (f"; failing: {', '.join(v['guard_failed'])}" if v["guard_failed"] else "") + ").", ""] + table(v["guard"])
    return "\n".join(L) + "\n"


def main() -> None:
    from .report import ROOT
    args = sys.argv[1:]
    if args[:1] == ["validate"]:
        v = validate()
        (ROOT / "reports" / "barrel_validation.md").write_text(validation_md(v))
        (ROOT / "reports" / "barrel_validation.json").write_text(json.dumps(v, indent=1))
        print(validation_md(v))
        return
    for pdb, res in barrel_counts({a.upper(): None for a in args}).items():
        print(pdb, json.dumps(res))


if __name__ == "__main__":
    main()
