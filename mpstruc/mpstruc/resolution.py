"""Parse mpstruc's free-text <resolution> field.

Grammar (applied to the decoded text):
  NUMBER [Å] [(qualifier)]          -> value = NUMBER
  NUMBER [Å] / NUMBER [Å] ...       -> value = min(NUMBERs), kind "multiple"
  text mentioning NMR / solid-state -> value None, method hint
  text mentioning electron crystallography / cryo-EM -> value None, method hint
  anything else ('-', 'HEK293 cells', ...) -> value None, kind "unparseable"
The raw string is always kept by the caller; the method hint is only a hint
(the real method comes from the PDB's exptl.method).
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from .text import to_text

_NUM = r"(\d+(?:\.\d+)?)\s*(?:Å|A)?"
_SINGLE = re.compile(rf"{_NUM}\s*(?:\((?P<q>[^)]*)\))?\Z")
_MULTI = re.compile(rf"{_NUM}(?:\s*/\s*{_NUM})+\s*(?:\((?P<q>[^)]*)\))?\Z")


@dataclass(frozen=True)
class Resolution:
    raw: str | None
    value: float | None
    kind: str  # "number" | "multiple" | "method" | "missing" | "unparseable"
    qualifier: str | None = None
    method_hint: str | None = None


def parse_resolution(raw: str | None) -> Resolution:
    text = to_text(raw)
    if text is None:
        return Resolution(raw, None, "missing")
    if m := _SINGLE.match(text):
        return Resolution(raw, float(m[1]), "number", m["q"])
    if m := _MULTI.match(text):
        nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", text.split("(")[0])]
        return Resolution(raw, min(nums), "multiple", m["q"])
    low = text.lower()
    if "nmr" in low or "strucuture" in low:
        hint = "solid-state NMR" if "solid" in low else "NMR"
        return Resolution(raw, None, "method", method_hint=hint)
    if "electron crystallography" in low:
        return Resolution(raw, None, "method", method_hint="electron crystallography")
    if "cryo" in low or re.search(r"\bem\b", low):
        return Resolution(raw, None, "method", method_hint="cryo-EM")
    return Resolution(raw, None, "unparseable")
