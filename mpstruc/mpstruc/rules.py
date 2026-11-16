"""The curators' rules for placing structures in a subfamily (curator_rules.json).

A rule names a subfamily as mpstruc does and quotes the curator's sentence:

    {"subgroup": "Outer Membrane Autotransporters", "subgroupKey": "G02-S006",
     "note": "Typically 12-stranded beta-barrel.", "barrelStrands": 12, "tolerance": 1,
     "addedBy": "curator", "addedOn": "2026-11-13"}

The subfamily is matched by its name in the mpstruc snapshot in use (keys shift
when subfamilies are added); the key in the file is kept as a fallback and for
readers. "barrelStrands" (optional) is the number of strands of the membrane
beta-barrel its members typically have, "tolerance" how far a counted number
may be from it and still fit (default 0). More rules of this kind are added by
editing the file; other criteria can be added as new optional fields, which
this loader passes through untouched.

automation.py uses the rules when it suggests a subfamily (apply_rules) and
copies them, with the official mpstruc entries that contradict one, to
../data/curator_rules.json for the site.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

RULES_PATH = Path(__file__).resolve().parent.parent / "curator_rules.json"
_REQUIRED = ("subgroup", "note")


class RuleError(ValueError):
    pass


def _norm(name: str) -> str:
    return re.sub(r"\s+", " ", name or "").strip().lower()


def load_rules(path: Path | None = None) -> list[dict]:
    """The rules in the file, checked. A missing file means no rules."""
    path = path or RULES_PATH
    if not path.exists():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    rules = data.get("rules") if isinstance(data, dict) else data
    if not isinstance(rules, list):
        raise RuleError(f"{path.name}: 'rules' must be a list")
    out = []
    for n, r in enumerate(rules, 1):
        if not isinstance(r, dict):
            raise RuleError(f"{path.name}: rule {n} must be an object")
        for k in _REQUIRED:
            if not isinstance(r.get(k), str) or not r[k].strip():
                raise RuleError(f"{path.name}: rule {n} needs a '{k}'")
        strands = r.get("barrelStrands")
        if strands is not None and (isinstance(strands, bool) or not isinstance(strands, int) or not 4 <= strands <= 60):
            raise RuleError(f"{path.name}: rule {n}: barrelStrands must be a whole number from 4 to 60")
        tol = r.get("tolerance", 0)
        if isinstance(tol, bool) or not isinstance(tol, int) or tol < 0:
            raise RuleError(f"{path.name}: rule {n}: tolerance must be 0 or a positive whole number")
        out.append({**r, "tolerance": tol})
    return out


def resolve(rules: list[dict], sg_names: dict[str, str]) -> dict[str, dict]:
    """{subgroup key in this snapshot: rule}. sg_names: {key: name as displayed}.

    Matched by name (case and spacing ignored); if no subfamily has that name, by the
    recorded key. A rule that matches nothing is left out (and reported by unmatched)."""
    by_name: dict[str, list[str]] = {}
    for key, name in sg_names.items():
        by_name.setdefault(_norm(name), []).append(key)
    out: dict[str, dict] = {}
    for r in rules:
        keys = by_name.get(_norm(r["subgroup"])) or []
        if len(keys) > 1 and r.get("subgroupKey") in keys:  # two subfamilies share a name: the recorded one
            keys = [r["subgroupKey"]]
        key = keys[0] if len(keys) == 1 else r.get("subgroupKey") if not keys and r.get("subgroupKey") in sg_names else None
        if key:
            out[key] = {**r, "key": key, "name": sg_names[key]}
    return out


def unmatched(rules: list[dict], sg_names: dict[str, str]) -> list[str]:
    """Subfamily names in the rules that this snapshot does not have."""
    got = {_norm(r["subgroup"]) for r in resolve(rules, sg_names).values()}
    return [r["subgroup"] for r in rules if _norm(r["subgroup"]) not in got]


def strand_rules(resolved: dict[str, dict]) -> dict[str, dict]:
    """The rules that give a number of barrel strands."""
    return {k: r for k, r in resolved.items() if r.get("barrelStrands")}


def fits(rule: dict, count: int | None) -> bool:
    """Does a counted number of barrel strands fit the rule?"""
    return count is not None and abs(count - rule["barrelStrands"]) <= rule.get("tolerance", 0)


def verdict(rule: dict, count: int, closed: bool = True, others=(), up_to: int | None = None) -> str:
    """"fits", "differs" or "unclear" for a barrel counted in a structure.

    A closed barrel fits or differs. If a larger ring fits the same strands (up_to, see
    barrel.ring_detours) and only one of the two numbers fits, it is "unclear".

    The count of a barrel whose seam is open (a barrel caught while it folds or with a
    substrate in it) is a lower bound: strands may be missing from it. More strands than
    the rule allows still contradict it ("differs"); fewer say nothing ("unclear"); and a
    number in the rule's range fits only if no other rule (others: the other strand rules)
    with more strands is still within reach, else an under-counted 16-strand Omp85 barrel
    would be confirmed as a 12-strand autotransporter ("unclear")."""
    lo, hi = rule["barrelStrands"] - rule.get("tolerance", 0), rule["barrelStrands"] + rule.get("tolerance", 0)
    if closed:
        if up_to and up_to > count and fits(rule, count) != fits(rule, up_to):
            return "unclear"
        return "fits" if fits(rule, count) else "differs"
    if count > hi:
        return "differs"
    if count < lo:
        return "unclear"
    reach = any(o["barrelStrands"] > rule["barrelStrands"] and count <= o["barrelStrands"] + o.get("tolerance", 0)
                and not fits(rule, o["barrelStrands"]) for o in others)
    return "unclear" if reach else "fits"


def possible(rule: dict, count: int, closed: bool = True, up_to: int | None = None) -> bool:
    """Could a barrel counted like this have the rule's number of strands?"""
    if closed:
        return fits(rule, count) or bool(up_to and fits(rule, up_to))
    return count <= rule["barrelStrands"] + rule.get("tolerance", 0)


def strands_phrase(rule: dict) -> str:
    """'typically 16-stranded'"""
    return f"typically {rule['barrelStrands']}-stranded"


def site_data(resolved: dict[str, dict], conflicts: list[dict] | None = None, generated: str = "",
              validation: dict | None = None) -> dict:
    """What the site shows: the rules by subfamily key, the official mpstruc entries that
    contradict one, and how well the strand counting did on known structures."""
    keep = ("subgroup", "note", "barrelStrands", "tolerance", "addedBy", "addedOn")
    out = {"generated": generated,
           "rules": [{"key": k, "name": r["name"], **{f: r[f] for f in keep if r.get(f) is not None}}
                     for k, r in sorted(resolved.items())],
           "conflicts": conflicts or []}
    if validation:
        out["counting"] = validation
    return out
