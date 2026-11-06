"""Add mpstruc-style <sub>/<sup> markup to plain text.

mpstruc curators write "Ca<sup>2+</sup>", "cytochrome b<sub>6</sub>f",
"GABA<sub>A</sub>", "5-HT<sub>2A</sub>" and so on. Rather than guessing,
this module learns how they write it from the database itself:

* `build_lexicon()` scans every name, description and bibliography title in
  the mpstruc snapshots, collects each marked-up word ("unit") such as
  "G<sub>s</sub>", and pairs it with its plain spelling ("Gs").  It then counts
  how often that plain spelling appears in mpstruc *without* markup and keeps
  only spellings curators mark up consistently (so a protein name that
  mpstruc usually leaves plain, e.g. "PGH2", is never tagged).
* Two small general rules cover forms the lexicon has not seen yet:
  metal-ion charges ("Sr2+") and simple salts/molecules ("SrCl2"), limited to
  the element symbols mpstruc itself marks up that way.

`format_text(s)` applies both.  It is idempotent, never touches HTML tags,
entities, URLs or PDB IDs, and never changes a word that already has markup.

Regenerate the lexicon after a new snapshot:  python -m mpstruc.markup build
Measure accuracy on held-out data:           python -m mpstruc.markup eval
"""
from __future__ import annotations

import html
import json
import re
import sys
import zlib
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "fixtures"
LEXICON_PATH = ROOT / "style_lexicon.json"
SNAPSHOTS = ("mpstruc_2026-10-02.xml", "mpstruc_2026-10-15.xml")

# Lexicon admission thresholds.  They were chosen by looking at the results
# of `evaluate`, so its figures are validation results, not an untouched test.
MIN_MARKED = 2          # seen marked-up at least this often
MIN_ENTRIES = 2         # ... in at least this many different entries
MIN_CONSISTENCY = 0.9   # share of all occurrences written in the chosen form
# Short keys such as H1, H4, B4, Y1, P63, ETA or c10 mean something else in
# other text (a helix, a histone, a residue, a ligand code).  A key of this
# length or shorter, learned from fewer than CONTEXT_FREE_ENTRIES entries, is
# only applied next to a word curators used it with ("histamine H4",
# "leukotriene B4", "space group P63", "ETA receptor").
CONTEXT_MAX_LEN = 3
CONTEXT_FREE_ENTRIES = 6

# ---------------------------------------------------------------- tokenising

_WORDCH = r"(?:[A-Za-z0-9Ͱ-Ͽ]|&[A-Za-z]+;)"
_TAGGED = r"<(sub|sup)>([^<]{1,12})</\1>"
# A marked-up unit: optional word characters, then one or more <sub>/<sup>
# runs each optionally followed by more word characters.
# A lipid-style locant prefix such as "PI(3,5)" belongs to the unit too.
UNIT_RE = re.compile(rf"(?:[A-Za-z]+\([0-9,]+\))?{_WORDCH}*(?:{_TAGGED}{_WORDCH}*)+")
SUBSUP_RE = re.compile(r"<(sub|sup)>(.*?)</\1>", re.S | re.I)
TAG_RE = re.compile(r"</?[A-Za-z][^>]*>")
URL_RE = re.compile(r"(?:https?://|www\.|doi:\s*|10\.\d{4,}/)\S+", re.I)
ENTITY_RE = re.compile(r"&(#[0-9]+|#[xX][0-9a-fA-F]+|[A-Za-z][A-Za-z0-9]*);")
# PDB IDs: a digit then three alphanumerics including at least one letter.
PDB_RE = re.compile(r"(?<![A-Za-z0-9])[0-9](?=[A-Za-z0-9]{0,2}[A-Za-z])[A-Za-z0-9]{3}(?![A-Za-z0-9])")
GREEK_ENTITIES = {n: html.unescape(f"&{n};") for n in (
    "alpha beta gamma delta epsilon zeta eta theta iota kappa lambda mu nu xi "
    "omicron pi rho sigma tau upsilon phi chi psi omega").split()}
GREEK_ENTITIES.update({n.capitalize(): v.upper() for n, v in list(GREEK_ENTITIES.items())})

_ISWORD = re.compile(r"[A-Za-z0-9Ͱ-Ͽ]")
BLOCK = "\x01"     # masked text that also counts as a word character
BOUNDARY = "\x00"  # masked text that counts as a word boundary


def _is_word(ch: str) -> bool:
    return bool(ch) and (ch == BLOCK or bool(_ISWORD.match(ch)))


def strip_subsup(s: str) -> str:
    """Remove <sub>/<sup> tags but keep their text (Ca<sup>2+</sup> -> Ca2+)."""
    return re.sub(r"</?(?:sub|sup)>", "", s, flags=re.I)


def _canon(s: str) -> str:
    """Greek entities to characters, so &beta;2 and β2 share a lexicon key."""
    return re.sub(r"&([A-Za-z]+);", lambda m: GREEK_ENTITIES.get(m[1], m[0]), s)


def units(raw: str):
    """Yield (plain, marked) for every marked-up unit in a raw mpstruc field."""
    for m in UNIT_RE.finditer(raw):
        start = m.start()
        if start and _is_word(raw[start - 1]):
            continue
        marked = m.group(0)
        if any(c.strip() != c for _, c in SUBSUP_RE.findall(marked)):
            continue                      # "G<sub>s </sub>": stray spaces
        if marked.startswith("<"):        # "<sup>175</sup>TLESWSM": no base
            continue
        yield _canon(strip_subsup(marked)), _canon(marked)


# ---------------------------------------------------------------- corpus

FIELDS = ("name", "description", "title")
CURATOR_FIELDS = ("name", "description")


def corpus(paths=None, exclude=()):
    """Yield (pdbCode, field, raw) for every text field, de-duplicated across
    snapshots.  `exclude` is a set of PDB codes to leave out."""
    from .model import load
    seen = set()
    for p in paths or [FIXTURES / s for s in SNAPSHOTS]:
        db = load(Path(p))
        for e in db.entries():
            if e.pdbCode in exclude:
                continue
            vals = [("name", e.name), ("description", e.description),
                    ("title", e.bibliography.title)]
            vals += [("title", b.title) for b in e.secondary]
            for f, v in vals:
                if v and (e.pdbCode, f, v) not in seen:
                    seen.add((e.pdbCode, f, v))
                    yield e.pdbCode, f, v


def _plain_counts(texts, keys):
    """How often each key appears as a whole word in tag-stripped text."""
    counts = Counter()
    if not keys:
        return counts
    rx = _key_regex(tuple(sorted(keys, key=lambda k: (-len(k), k))))
    for t in texts:
        plain = _canon(strip_subsup(t))
        for m in rx.finditer(plain):
            if _bounded(plain, m.start(), m.end()):
                counts[m.group(0)] += 1
    return counts


def build_lexicon(records=None) -> dict:
    """Mine the corpus; return {"lexicon": {plain: marked}, "stats": ...}.

    Marked-up spellings are collected from every field.  Consistency is
    judged on the text curators write themselves (entry names and
    descriptions); bibliography titles are copied from the journals, which
    often leave "K+" or "A2A" plain, so they only decide for spellings that
    never occur in curator-written text."""
    records = list(records if records is not None else corpus())
    forms = defaultdict(Counter)
    where = defaultdict(set)
    for code, _, raw in records:
        for plain, marked in units(raw):
            forms[plain][marked] += 1
            where[plain].add(code)
    own = [r for r in records if r[1] in CURATOR_FIELDS]
    own_forms = defaultdict(Counter)
    for _, _, raw in own:
        for plain, marked in units(raw):
            own_forms[plain][marked] += 1
    own_totals = _plain_counts([r for _, _, r in own], forms.keys())
    all_totals = _plain_counts([r for _, _, r in records], forms.keys())
    lexicon, stats = {}, {}
    for plain, fc in forms.items():
        best, n_all = fc.most_common(1)[0]
        if own_totals[plain]:
            src, fcs, total = "curator", own_forms[plain], own_totals[plain]
        else:
            src, fcs, total = "all", fc, all_totals[plain]
        n = fcs.get(best, 0)
        total = max(total, sum(fcs.values()))
        consistency = n / total if total else 0.0
        keep = (n_all >= MIN_MARKED and len(where[plain]) >= MIN_ENTRIES
                and consistency >= MIN_CONSISTENCY and _plausible(plain, best))
        stats[plain] = {"form": best, "marked": n_all, "entries": len(where[plain]), "judged_on": src,
                        "marked_there": n, "seen_there": total,
                        "consistency": round(consistency, 3), "kept": keep}
        if keep:
            lexicon[plain] = best
    context = _learn_context(records, {k: v for k, v in lexicon.items()
                                       if _needs_context(k, stats[k]["entries"])})
    return {"lexicon": dict(sorted(lexicon.items())), "context": context,
            "stats": dict(sorted(stats.items()))}


def _needs_context(plain: str, entries: int) -> bool:
    return (len(plain) <= CONTEXT_MAX_LEN and entries < CONTEXT_FREE_ENTRIES
            and not re.search(r"[+\-−]", plain))


_PREV_WORD = re.compile(r"([A-Za-zͰ-Ͽ][A-Za-z0-9Ͱ-Ͽ]*)[ \t\n/(-]{1,3}$")
_NEXT_WORD = re.compile(r"[ \t\n/)-]{1,3}([A-Za-zͰ-Ͽ][A-Za-z0-9Ͱ-Ͽ]*)")


# Too common to tell one meaning of a short key from another.
_STOPWORDS = {"a", "an", "the", "and", "or", "with", "in", "of", "on", "at", "to",
              "by", "for", "from", "as", "is", "bound", "without"}


def _ctx_word(w: str) -> str:
    w = w.lower()
    w = w[:-1] if len(w) > 3 and w.endswith("s") else w
    return "" if w in _STOPWORDS else w


def _neighbours(s: str, i: int, j: int):
    """(word before, word after) a match, normalised; '' if none."""
    a = _PREV_WORD.search(s[max(0, i - 40):i])
    b = _NEXT_WORD.match(s, j)
    return (_ctx_word(a[1]) if a else "", _ctx_word(b[1]) if b else "")


def _learn_context(records, keys: dict) -> dict:
    """For each short key, the words written right before and after it in
    the fields where curators marked it up."""
    ctx = {k: {"before": set(), "after": set()} for k in keys}
    if not keys:
        return {}
    rx = _key_regex(tuple(sorted(keys, key=lambda k: (-len(k), k))))
    for _, _, raw in records:
        canon = _canon(raw)
        plain = _canon(strip_subsup(raw))
        for m in rx.finditer(plain):
            k = m.group(0)
            if keys[k] not in canon or not _bounded(plain, m.start(), m.end()):
                continue
            before, after = _neighbours(plain, m.start(), m.end())
            if before:
                ctx[k]["before"].add(before)
            if after:
                ctx[k]["after"].add(after)
    return {k: {side: sorted(v[side]) for side in v} for k, v in sorted(ctx.items())}


# Lexicon keys that are also English words: tagged only mid-sentence.
ENGLISH_WORDS = {"Go", "Be", "No", "So", "Do", "As", "In", "Is", "It", "On", "Or", "To", "By", "Me"}


def _plausible(plain: str, marked: str) -> bool:
    """Reject units too generic to replay safely in new text."""
    if len(plain) < 2 or not re.search(r"[A-Za-zͰ-Ͽ]", plain):
        return False
    if PDB_RE.fullmatch(plain):
        return False
    return True


# ---------------------------------------------------------------- general rules

# Element symbols mpstruc writes with a superscript charge.  A hand-written
# list (chosen while looking at the whole corpus), not mined by `build`.
CATIONS = ("Ag Al Ba Ca Cd Co Cs Cu Fe Gd Hg K Li Mg Mn Na Ni Pb Rb Sr Tb Tl Zn").split()
# "Sr2+", "Co3+", "Cs+": cation with charge.
CATION_RE = re.compile(rf"({'|'.join(sorted(CATIONS, key=len, reverse=True))})([1-4]?\+)")

# The simple-formula rule only accepts recognised chemistry, never "any
# string of element letters with a number": gene and protein names such as
# PHO84, OPN4, OSH6, SOP4, SCO1, NOS1, COS7 or ligand codes such as CO89 are
# made of element letters too.  Accepted are
#   * a known anion or small molecule on its own: SO4, HCO3-, AlF4-, H2O, CO2;
#   * a salt: metal (or ammonium) with an optional count of 2-4, then a
#     halide or one of the anions with an optional count of 2-4: MgCl2,
#     KH2PO4, Na2SO4, NH4Cl, AlCl3.
ANIONS = ("SO4 SO3 PO4 HPO4 H2PO4 NO3 CO3 HCO3 C2O4 ClO4 VO4 AlF4 BeF3 "
          "MgF3 MgF4").split()
MOLECULES = ("H2O H2O2 D2O CO2 PCO2 N2O H2S NH3 NH4 Mn4CaO5").split()
HALIDES = ("Cl", "Br", "F")
_alt = lambda xs: "|".join(sorted(xs, key=len, reverse=True))
SALT_RE = re.compile(rf"(?:{_alt(CATIONS + ['Be', 'NH4'])})(?:[2-4])?"
                     rf"(?:{_alt(HALIDES + tuple(ANIONS))})(?:[2-4])?")
# Formulas that look chemical but are names or codes mpstruc writes plain.
FORMULA_BLOCK = {"CO2+", "NO", "CN", "HCN", "NO2"}


def _formula_markup(tok: str) -> str | None:
    """'MgCl2' -> 'MgCl<sub>2</sub>', 'NH4+' -> 'NH<sub>4</sub><sup>+</sup>'.
    Returns None for anything that is not a recognised formula (see above)."""
    if tok in FORMULA_BLOCK or not re.fullmatch(r"[A-Za-z0-9]+[+-]?", tok):
        return None
    # "AlF4-": the 4 is a count.  "SO42-", "MgF42-", "C2O42-": a known anion
    # followed by 2-4 and a sign -- the last digit is the charge.
    charge = None
    if tok[-1] in "+-":
        m = re.fullmatch(r"(.*?)([2-4]?[+-])", tok)
        if m and (m[1] in ANIONS or m[1] == "NH4"):
            charge = m[2]
            body = m[1]
        elif tok[:-1] in ANIONS or tok[:-1] == "NH4":
            charge, body = tok[-1], tok[:-1]
        else:
            return None
        if charge[0].isdigit() and body == "NH4":
            return None
    else:
        body = tok
    if not re.search(r"\d", body):
        return None
    if not (body in ANIONS or body in MOLECULES or (not charge and SALT_RE.fullmatch(body))):
        return None
    out = re.sub(r"(?<=[A-Za-z])([0-9]+)", r"<sub>\1</sub>", body)
    return out + (f"<sup>{charge}</sup>" if charge else "")


# ---------------------------------------------------------------- applying

_LB = r"(?<![A-Za-z0-9\u0370-\u03ff\x01&])"
_RB = r"(?![A-Za-z0-9\u0370-\u03ff\x01])"


def _key_regex(keys: tuple[str, ...]):
    """One regex for all keys (longest first), with word boundaries built in
    so a failed boundary falls back to a shorter key instead of skipping."""
    return re.compile(_LB + "(?:" + "|".join(re.escape(k) for k in keys) + ")" + _RB)


def _bounded(s: str, i: int, j: int) -> bool:
    """Whole-unit match: no word character (or masked text) on either side,
    and a trailing '+'/'-' charge must not run on into another word."""
    if i > 0 and (_is_word(s[i - 1]) or s[i - 1] == "&"):
        return False
    if j < len(s) and (_is_word(s[j]) or s[j] == "+"):
        return False
    return True


@lru_cache(maxsize=1)
def _load():
    try:
        data = json.loads(LEXICON_PATH.read_text(encoding="utf-8"))
        data["lexicon"], data["context"]
    except (OSError, ValueError, KeyError):
        data = build_lexicon()
    return _prepare(data)


def _prepare(lex: dict):
    """Add entity/character variants of Greek keys; compile one regex.
    `lex` is either a plain {plain: marked} dict or the output of
    `build_lexicon` (which also carries the context words for short keys).
    Returns (keys -> markup, regex, keys -> context)."""
    context = {}
    if isinstance(lex.get("lexicon"), dict):
        context = lex.get("context") or {}
        lex = lex["lexicon"]
    full = dict(lex)
    inv = {v: k for k, v in GREEK_ENTITIES.items() if k.islower() or k[0].isupper()}
    for plain, marked in lex.items():
        if re.search("[Ͱ-Ͽ]", plain):
            ent = lambda s: re.sub("[Ͱ-Ͽ]", lambda m: f"&{inv[m[0]]};" if m[0] in inv else m[0], s)
            full.setdefault(ent(plain), ent(marked))
    keys = tuple(sorted(full, key=lambda k: (-len(k), k)))
    ctx = {k: context[_canon(k)] for k in full if _canon(k) in context}
    return full, _key_regex(keys) if keys else None, ctx


def _mask(s: str) -> str:
    """Same-length copy of s with everything we must not touch hidden.
    Existing <sub>/<sup> spans become word-like (so neighbours are left
    alone); other tags and non-Greek entities become boundaries."""
    out = list(s)

    def fill(a, b, ch):
        out[a:b] = ch * (b - a)

    for m in SUBSUP_RE.finditer(s):
        fill(m.start(), m.end(), BLOCK)
    for m in TAG_RE.finditer(s):
        if out[m.start()] != BLOCK:
            fill(m.start(), m.end(), BOUNDARY)
    for m in URL_RE.finditer(s):
        fill(m.start(), m.end(), BLOCK)
    for m in ENTITY_RE.finditer(s):
        if m[1] not in GREEK_ENTITIES and out[m.start()] == "&":
            fill(m.start(), m.end(), BOUNDARY)
    masked = "".join(out)
    for m in PDB_RE.finditer(masked):
        # "2Ca2+", "2NO3−" (two ions) look like PDB IDs; a real ID is never
        # followed by a charge sign.
        nxt = masked[m.end():m.end() + 2]
        if nxt[:1] in ("+", "−") or (nxt[:1] == "-" and not nxt[1:2].isalnum()):
            continue
        fill(m.start(), m.end(), BLOCK)
    return "".join(out)


_TOKEN_RE = re.compile(r"[A-Za-z0-9]+[+-]?")
# A chain of three-letter residue/ligand codes: "DAL-GLY-MEA-".
_CODE_CHAIN = re.compile(r"(?<![A-Za-z0-9])(?:[A-Z]{3}-){2,}$")


def format_text(s: str | None, lexicon: dict | None = None) -> str | None:
    """Return s with mpstruc-style <sub>/<sup> markup added.

    Plain words are matched against the mined lexicon first, then the
    cation-charge and simple-formula rules.  Anything already marked up, and
    any HTML tag, entity, URL or PDB ID, is left exactly as it was."""
    if not s:
        return s
    full, rx, ctx = _prepare(lexicon) if lexicon is not None else _load()
    masked = _mask(s)
    edits = []  # (start, end, replacement)
    taken = [False] * len(s)
    if rx is not None:
        for m in rx.finditer(masked):
            i, j = m.span()
            if not _bounded(masked, i, j):
                continue
            key = m.group(0)
            mk = full[key]
            if _CODE_CHAIN.search(masked[max(0, i - 40):i]):
                continue  # "TYR-DAL-GLY-MEA-ETA": a ligand code, not ETA
            pair = re.match(r"/([a-z]{1,5}|[0-9]{1,3}|[α-ω])(?![A-Za-z0-9Ͱ-Ͽ\x01])", masked[j:])
            if not (pair and re.search(r"G<sub>[^<]+</sub>$", mk)):
                pair = None
            if key in ctx:
                end = j + (len(pair.group(0)) if pair else 0)
                before, after = _neighbours(masked, i, end)
                if before not in ctx[key]["before"] and after not in ctx[key]["after"]:
                    continue  # "helix H1", "Histone H4": not the learned context
            if (key in ENGLISH_WORDS and re.search(r"(?:^|[.!?:]\s*|\x00\s*)$", masked[:i])
                    and not re.match(r"[-/]|\s+(?:protein|complex|heterotrimer)\b", masked[j:])):
                continue  # "Go" opening a sentence is the verb, not G-alpha-o
            if pair:
                # "Gs/q", "Gq/11", "G12/13": mpstruc subscripts the whole
                # pair (miniG<sub>s/q</sub>).
                mk = mk[:-len("</sub>")] + pair.group(0) + "</sub>"
                j += len(pair.group(0))
            edits.append((i, j, mk))
            taken[i:j] = [True] * (j - i)
    for m in _TOKEN_RE.finditer(masked):
        i, j = m.span()
        if any(taken[i:j]) or not _bounded(masked, i, j):
            continue
        tok = m.group(0)
        if _CODE_CHAIN.search(masked[max(0, i - 40):i]):
            continue
        c = CATION_RE.fullmatch(tok)
        if c:
            rep = f"{c[1]}<sup>{c[2]}</sup>"
        else:
            rep = _formula_markup(tok)
        if tok.endswith("-") and rep and j < len(s) and s[j:j + 1].isalnum():
            rep = None  # "Cl-dependent": a hyphen, not a charge
        if rep:
            edits.append((i, j, rep))
    for i, j, rep in sorted(edits, reverse=True):
        s = s[:i] + rep + s[j:]
    return s


# ---------------------------------------------------------------- measuring

HELD_OUT_19 = ("11CZ 11VV 13CH 22FR 3SQ6 8FYN 8HK3 8HK5 8QMP 8ZRP 9ETZ 9JEL "
               "9L0B 9WPQ 9WQA 9XN2 9XN3 9XN4 9Y1N").split()


def _tag_spans(raw: str) -> set:
    """Positions of sub/sup tags in plain-text coordinates: {(start, end, kind)}."""
    spans, pos, i = set(), 0, 0
    for m in SUBSUP_RE.finditer(raw):
        pos += len(raw[i:m.start()])
        n = len(m.group(2))
        spans.add((pos, pos + n, m.group(1).lower()))
        pos += n
        i = m.end()
    return spans


def score(records, lexicon) -> dict:
    """Strip markup from each field, re-add it, compare with mpstruc."""
    tp = fp = fn = near = 0
    changed_plain = exact = marked_fields = plain_fields = 0
    errors = []
    for code, f, raw in records:
        plain = strip_subsup(raw)
        out = format_text(plain, lexicon)
        gold, pred = _tag_spans(raw), _tag_spans(out)
        tp += len(gold & pred)
        fp += len(pred - gold)
        fn += len(gold - pred)
        # Wrong tags that still land on a word mpstruc marks up, just split
        # differently (Na<sub>v</sub> vs Na<sub>V</sub>, Ca2<sup>+</sup>).
        words = [m.span() for m in re.finditer(r"[^\s,;()]+", plain)]
        gold_words = {w for w in words for a, b, _ in gold if w[0] <= a < w[1]}
        near += sum(1 for a, b, _ in pred - gold
                    if any(w[0] <= a < w[1] for w in gold_words))
        if gold:
            marked_fields += 1
            exact += out == raw
        else:
            plain_fields += 1
            changed_plain += out != plain
        if pred - gold and len(errors) < 40:
            errors.append((code, f, [plain[a:b] for a, b, _ in sorted(pred - gold)],
                           re.sub(r"\s+", " ", out)[:160]))
    return {"tags_gold": tp + fn, "tags_added": tp + fp, "correct": tp,
            "precision": round(tp / (tp + fp), 4) if tp + fp else 1.0,
            "recall": round(tp / (tp + fn), 4) if tp + fn else 1.0,
            "wrong_split_on_marked_word": near,
            "marked_fields": marked_fields, "marked_fields_exact": exact,
            "plain_fields": plain_fields, "plain_fields_changed": changed_plain,
            "false_positive_examples": errors}


def _families(paths=None) -> dict:
    """PDB code -> PDB code of its top-level (parent) entry, so related
    structures from one paper land in the same fold."""
    from .model import load
    fam = {}
    for p in paths or [FIXTURES / s for s in SNAPSHOTS]:
        for e in load(Path(p)).entries():
            for d in [e] + list(e.daughters):
                fam.setdefault(d.pdbCode, e.pdbCode if e.kind == "parent" else fam.get(e.pdbCode, e.pdbCode))
    return fam


def _split_scores(test, lexicon) -> dict:
    return {"curator_fields": score([r for r in test if r[1] in CURATOR_FIELDS], lexicon),
            "bibliography_titles": score([r for r in test if r[1] not in CURATOR_FIELDS], lexicon)}


def evaluate(folds: int = 5) -> dict:
    """Two held-out tests, each reported separately for curator-written text
    (entry names and descriptions) and for bibliography titles (copied from
    the journals, so they follow the journal's style, not mpstruc's):
    1. k-fold cross-validation over the newest snapshot, split by entry
       family so related structures never sit on both sides;
    2. the 19 entries curators added between the two snapshots, with the
       lexicon learned from the older snapshot only.
    Both were used while tuning, so treat the figures as validation results."""
    recs = list(corpus())
    fam = _families()
    keys = ("tags_gold", "tags_added", "correct", "wrong_split_on_marked_word", "marked_fields",
            "marked_fields_exact", "plain_fields", "plain_fields_changed")
    agg = defaultdict(Counter)
    examples = defaultdict(list)
    for k in range(folds):
        fold = lambda c: zlib.crc32(fam.get(c, c).encode()) % folds == k
        train = [r for r in recs if not fold(r[0])]
        test = [r for r in recs if fold(r[0])]
        for part, res in _split_scores(test, build_lexicon(train)).items():
            for key in keys:
                agg[part][key] += res[key]
            examples[part] += res["false_positive_examples"]
    cv = {}
    for part, a in agg.items():
        cv[part] = dict(a)
        cv[part]["precision"] = round(a["correct"] / a["tags_added"], 4) if a["tags_added"] else 1.0
        cv[part]["recall"] = round(a["correct"] / a["tags_gold"], 4) if a["tags_gold"] else 1.0
        cv[part]["false_positive_examples"] = examples[part][:25]
    old = list(corpus([FIXTURES / SNAPSHOTS[0]]))
    lex_old = build_lexicon(old)
    new19 = [r for r in corpus([FIXTURES / SNAPSHOTS[1]]) if r[0] in HELD_OUT_19]
    return {"note": "Validation figures: the thresholds and the hand-written word lists "
                    "(CATIONS, ANIONS, MOLECULES, FORMULA_BLOCK, ENGLISH_WORDS) were chosen "
                    "while looking at these same results, so they are not an untouched test. "
                    "The next mpstruc snapshot is the honest test.",
            "cross_validation": cv, "new_entries_19": _split_scores(new19, lex_old)}


def build(path: Path = LEXICON_PATH) -> dict:
    data = build_lexicon()
    data = {"about": "Plain spelling -> mpstruc markup, mined from mpstruc "
                     "names, descriptions and bibliography titles "
                     f"({', '.join(SNAPSHOTS)}). Regenerate: python -m mpstruc.markup build",
            "min_marked": MIN_MARKED, "min_consistency": MIN_CONSISTENCY, **data}
    path.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    _load.cache_clear()
    return data


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    cmd = argv[0] if argv else "help"
    if cmd == "build":
        d = build()
        print(f"wrote {LEXICON_PATH} ({len(d['lexicon'])} forms kept of {len(d['stats'])} seen)")
    elif cmd == "eval":
        print(json.dumps(evaluate(), indent=1, ensure_ascii=False))
    elif cmd == "format":
        print(format_text(" ".join(argv[1:])))
    else:
        print("usage: python -m mpstruc.markup build | eval | format TEXT")


if __name__ == "__main__":
    main()
