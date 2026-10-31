"""Turn raw mpstruc field markup into display text.

Raw fields mix HTML entities, a few known-bad entities, <br />, inline tags
(<sub>, <sup>, <em>, <i>, <b>, <a>) and some malformed tags. These helpers
never modify stored raw text; they only derive views of it.
"""
from __future__ import annotations

import html
import re
from html.entities import html5

# Misspelled entity names seen in the source, and what the curator meant.
ENTITY_REPAIRS = {"Ariing": "Å", "alfa": "α"}
# Broken spellings of &Aring; that are not entity references at all:
# "& Aring;" (space after '&') and "Aring;" (no '&').
BROKEN_ARING_RE = re.compile(r"&\s+Aring;|(?<![&\w])Aring;")
# Invisible format characters (zero-width space/joiners, word joiner, BOM).
INVISIBLE_RE = re.compile("[\u200b\u200c\u200d\u2060\ufeff]")
# "&uml;" inside a word, e.g. "M&uml;ller": a diaeresis typed instead of &uuml;.
STRAY_UML_RE = re.compile(r"(?<=[A-Za-z])&uml;(?=[A-Za-z])")

ENTITY_RE = re.compile(r"&(#[0-9]+|#[xX][0-9a-fA-F]+|[A-Za-z][A-Za-z0-9]*);")
# An entity name not followed by ';' (e.g. "&Aring:" or "&Ouml,").
UNTERMINATED_RE = re.compile(r"&([A-Za-z][A-Za-z0-9]*)(?![A-Za-z0-9;])")
TAG_RE = re.compile(r"</?[A-Za-z][^>]*>")
BR_RE = re.compile(r"<br\s*/?>", re.I)
PDB_ID_RE = re.compile(r"\b[0-9][A-Za-z0-9]{3}\b")


def decode_entities(s: str) -> str:
    """Decode as a browser would (html.unescape, which also accepts legacy
    entities with no ';'), after repairing the known misspellings."""
    for bad, good in ENTITY_REPAIRS.items():
        s = s.replace(f"&{bad};", good)
    s = BROKEN_ARING_RE.sub("Å", s)
    return html.unescape(s)


def unknown_entities(s: str | None) -> list[str]:
    """Entity names that are neither HTML5 nor numeric (e.g. Ariing, alfa)."""
    if not s:
        return []
    return [m[1] for m in ENTITY_RE.finditer(s)
            if not m[1].startswith("#") and m[1] + ";" not in html5]


def broken_entities(s: str | None) -> list[str]:
    """Malformed entity spellings repaired or left visible in display text."""
    if not s:
        return []
    return [m[0] for m in BROKEN_ARING_RE.finditer(s)] + [m[0] for m in STRAY_UML_RE.finditer(s)]


def unterminated_entities(s: str | None) -> list[str]:
    """Entity references missing their ';' (browsers guess; XML rejects)."""
    if not s:
        return []
    return [m[0] for m in UNTERMINATED_RE.finditer(s) if m[1] + ";" in html5]


def to_text(raw: str | None) -> str | None:
    """Plain text: <br/> becomes a newline, other tags are dropped, entities decoded."""
    if raw is None:
        return None
    s = BR_RE.sub("\n", raw)
    s = TAG_RE.sub("", s)
    s = INVISIBLE_RE.sub("", decode_entities(s).replace("\xa0", " "))
    lines = [re.sub(r"[ \t]+", " ", ln).strip() for ln in s.split("\n")]
    s = "\n".join(lines).strip("\n")  # keep paragraph breaks (<br /><br />)
    return s or None


def mentioned_pdb_ids(raw: str | None) -> list[str]:
    """PDB IDs cited inline in mpstruc's style, right after a resolution:
    '2.22 Å: 6BL3', '1.85 Å 6S63', '2.7 Å, 1KF6', '2.36 Å; 8E5V', '3.30 Å7V2D',
    '2.15 Å: : 1PO0', or a list such as 'Å: 6BL3, 6BL4 & 6BL5'. Requiring
    the Å context avoids matching concentrations ('10mM'), mutations ('296A')
    or notes about former codes. IDs are upper-cased; callers should check
    them against the PDB.
    """
    text = to_text(raw) or ""
    out: list[str] = []
    for m in re.finditer(r"Å[\s:;,]*((?:[0-9][A-Za-z0-9]{3}\b(?:\s*(?:,|&|and)\s*)?)+)", text):
        for pid in PDB_ID_RE.findall(m[1]):
            pid = pid.upper()
            if pid not in out:
                out.append(pid)
    return out
