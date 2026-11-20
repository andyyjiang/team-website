"""Export the database as one JSON file for the static web page.

    python -m mpstruc.web [--offline]
    python -m mpstruc.web links [--mpstruc FILE] [--candidates FILE]

Writes ../data/mpstruc.json (relative to this project folder), which the
site's mpstruc.html page loads. `links` only refreshes the EMDB map IDs
("emdb") of every structure in the existing ../data/mpstruc.json and
../data/candidates.json, in place (see patch_links). Curated fields are kept as sanitised HTML
(sub/sup/italics survive; broken tags are dropped); PDB-derived fields are
kept separate so the page can label them "from PDB".
"""
from __future__ import annotations

import argparse
import html
import json
import re
from pathlib import Path

from . import integrity
from .model import Database, Entry, load
from .report import source_xml
from .resolution import parse_resolution
from .text import BR_RE, INVISIBLE_RE, decode_entities, to_text

OUT = Path(__file__).resolve().parent.parent.parent / "data" / "mpstruc.json"
CANDIDATES = OUT.parent / "candidates.json"

ALLOWED = {"sub", "sup", "em", "i", "b"}
RENAME = {"super": "sup"}
PDB_LINK_RE = re.compile(r"(Å[\s:;,]*)((?:[0-9][A-Za-z0-9]{3}\b(?:\s*(?:,|&amp;|and)\s*)?)+)")


def safe_html(raw: str | None, link_pdb: bool = False) -> str | None:
    """Raw mpstruc markup -> HTML with only <br>, <sub>, <sup>, <em>, <i>, <b>.

    Text is entity-decoded and re-escaped; unknown or malformed tags
    (<ahref=...>, <a ...>, <super> is renamed) are dropped, and unclosed
    tags are closed. With link_pdb, PDB IDs cited after 'Å' become links.
    """
    if raw is None:
        return None
    out: list[str] = []
    stack: list[str] = []
    for part in re.split(r"(<[^>]*>)", BR_RE.sub("\n", raw)):
        if not part:
            continue
        if part.startswith("<"):
            m = re.match(r"<\s*(/?)\s*([A-Za-z]+)", part)
            if not m:
                continue
            closing, name = m[1] == "/", RENAME.get(m[2].lower(), m[2].lower())
            if name not in ALLOWED:
                continue
            if closing:
                if name in stack:
                    while stack:
                        t = stack.pop()
                        out.append(f"</{t}>")
                        if t == name:
                            break
            else:
                stack.append(name)
                out.append(f"<{name}>")
            continue
        out.append(html.escape(INVISIBLE_RE.sub("", decode_entities(part)).replace("\xa0", " "), quote=False))
    out += [f"</{t}>" for t in reversed(stack)]
    s = "".join(out)
    lines = [re.sub(r"[ \t]+", " ", ln).strip() for ln in s.split("\n")]
    s = "<br>".join(lines).strip()
    s = re.sub(r"^(<br>)+|(<br>)+$", "", s)
    s = re.sub(r"<(sub|sup|em|i|b)></\1>", "", s)
    if link_pdb:  # after joining, so a dropped tag between 'Å:' and the ID doesn't block the link
        s = PDB_LINK_RE.sub(lambda m: m[1] + re.sub(
            r"\b([0-9][A-Za-z0-9]{3})\b",
            lambda n: f'<a class="pdb-link" href="#/pdb/{n[1].upper()}">{n[1]}</a>', m[2]), s)
    return s or None


# Expression hosts: (pattern, canonical name). First match wins.
HOSTS = [
    (r"cell[- ]?free", "Cell-free expression"),
    (r"\be\.?\s*coli\b|escherichia", "Escherichia coli"),
    (r"\bhek|human embryonic kidney|expi293|\b293\b|293[ft-]", "HEK293 cells"),
    (r"spodoptera|\bsf[- ]?9\b|\bsf[- ]?21\b|s\.\s*frugiperda|insect cells|baculovirus", "Insect cells (Spodoptera frugiperda)"),
    (r"trichoplusia|high[- ]?five|\bhi5\b|t\.\s*ni\b", "Insect cells (Trichoplusia ni)"),
    (r"pichia|koma|kmoag|k\.\s*pastoris|p\.\s*pastoris|k\.\s*phaffii", "Pichia pastoris (Komagataella)"),
    (r"saccharomyces|s\.\s*cerevisiae|yeast", "Saccharomyces cerevisiae"),
    (r"\bcho\b|chinese hamster", "CHO cells"),
    (r"lactococcus|l\.\s*lactis", "Lactococcus lactis"),
    (r"xenopus|oocyte", "Xenopus oocytes"),
    (r"drosophila|\bs2\b", "Drosophila S2 cells"),
    (r"smegmatis", "Mycobacterium smegmatis"),
    (r"bacillus subtilis|b\.\s*subtilis", "Bacillus subtilis"),
    (r"pseudomonas", "Pseudomonas"),
    (r"leishmania|lexsy", "Leishmania tarentolae"),
    (r"hela", "HeLa cells"),
    (r"cos[- ]?\d|cos cells", "COS cells"),
    (r"bhk", "BHK cells"),
    (r"synthes|chemical|peptide synth", "Chemical synthesis"),
]


def norm_host(raw: str | None) -> str:
    t = (to_text(raw) or "").strip()
    if not t:
        return "Not stated (native source)"
    low = t.lower()
    for pat, name in HOSTS:
        if re.search(pat, low):
            return name
    return t


def _bib(b) -> dict:
    return {k: v for k, v in {
        "authors": to_text(b.authors), "year": to_text(b.year), "title": safe_html(b.title),
        "journal": to_text(b.journal), "volume": to_text(b.volume), "issue": to_text(b.issue),
        "pages": to_text(b.pages), "doi": to_text(b.doi), "pmid": to_text(b.pubMedId),
        "notes": safe_html(b.notes),
    }.items() if v}


def _method(methods: list[str]) -> str | None:
    if not methods:
        return None
    m = " + ".join(sorted({x.title().replace("Nmr", "NMR").replace("X-Ray", "X-ray") for x in methods}))
    return m.replace("Electron Microscopy", "Cryo-EM")


def _pdb_info(s: dict | None, emdb: list[str] | None = None) -> dict | None:
    """PDB facts shown on the page; "emdb" lists the entry's EMDB maps (only when it has any)."""
    if not s:
        return None
    res = (s.get("rcsb_entry_info") or {}).get("resolution_combined") or []
    return {k: v for k, v in {
        "title": (s.get("struct") or {}).get("title"),
        "method": _method([m.get("method", "") for m in s.get("exptl") or []]),
        "res": min(res) if res else None,
        "released": ((s.get("rcsb_accession_info") or {}).get("initial_release_date") or "")[:10] or None,
        "emdb": list(emdb) if emdb else None,
    }.items() if v is not None}


def _entry(e: Entry, flags: dict, pdb: dict, emdb: dict | None = None) -> dict:
    emdb = emdb or {}
    res = parse_resolution(e.resolution)
    d = {
        "k": e.key, "pdb": e.pdbCode, "kind": e.kind[0],  # "p" | "d"
        "name": safe_html(e.name), "species": safe_html(e.species),
        "domain": to_text(e.taxonomicDomain), "host": to_text(e.expressedInSpecies),
        "hostN": norm_host(e.expressedInSpecies), "res": to_text(e.resolution), "resV": res.value,
        "resHint": res.method_hint, "desc": safe_html(e.description, link_pdb=True),
        "year": int(y) if (y := to_text(e.bibliography.year)) and y.isdigit() else None,
        "bib": [_bib(b) for b in (e.bibliography, *e.secondary)],
        "rel": [r.pdb for r in e.related],
        "x": _pdb_info(pdb.get(e.pdbCode), emdb.get(e.pdbCode)),
        "flags": flags.get(e.key),
    }
    if e.kind == "parent":
        d["d"] = [_entry(x, flags, pdb, emdb) for x in e.daughters]
    return {k: v for k, v in d.items() if v not in (None, [], "")}


def _counts(db: Database) -> dict:
    ents = list(db.entries())
    parents = {e.pdbCode.upper() for e in ents if e.kind == "parent"}
    members = {e.pdbCode.upper() for e in ents if e.kind == "daughter"} - parents
    related = {r.pdb.upper() for e in ents for r in e.related} - parents - members
    return {
        "parents": len(parents), "daughters": len(members | related),
        "daughtersListed": len(members), "daughtersRelated": len(related),
        "structures": len(parents | members | related),
        "listings": len(ents) + sum(len(e.related) for e in ents),
        "subgroups": len(list(db.subgroups())),
    }


def build(db: Database, online: bool = True) -> dict:
    findings = integrity.run(db, online=online)
    flags: dict[str, list] = {}
    for f in findings:
        if f["severity"] == "info" or not f["key"]:
            continue
        owner = f["key"] if "-R" not in f["key"] else f["key"].rsplit("-R", 1)[0]
        msg = f["message"] if owner == f["key"] else f"Related {f['pdb']}: {f['message']}"
        flags.setdefault(owner, []).append([f["severity"][0], f["check"], msg])

    ents = list(db.entries())
    related = sorted({r.pdb for e in ents for r in e.related})
    pdb: dict = {}
    emdb: dict = {}
    rel_info: dict = {}
    if online:
        from . import rcsb
        pdb = rcsb.entry_summaries([e.pdbCode for e in ents])
        try:  # optional: without it the page still links to the PDB; the "links" step fills the maps in later
            emdb = rcsb.emdb_ids(sorted({e.pdbCode for e in ents} | set(related)))
        except Exception as ex:  # noqa: BLE001
            print(f"::warning::EMDB maps not fetched ({type(ex).__name__}: {ex}); building without them")
            emdb = {}
        for pid, s in rcsb.entry_summaries(related).items():
            if info := _pdb_info(s, emdb.get(pid)):
                rel_info[pid] = info

    groups = []
    for g in db.groups:
        sgs = []
        for sg in g.subgroups:
            parents = [_entry(p, flags, pdb, emdb) for p in sg.parents]
            n = sum(1 + len(p.daughters) + len(p.related) + sum(len(d.related) for d in p.daughters)
                    for p in sg.parents)
            sgs.append({"k": sg.key, "name": safe_html(sg.name), "n": n, "parents": parents})
        groups.append({"k": g.key, "name": g.name, "topology": g.topology, "subgroups": sgs})

    meta = db.meta
    return {
        "meta": {
            "timestamp": meta.get("timeStamp"), "lastNewProtein": meta.get("lastNewProteinDate"),
            "lastEdit": meta.get("lastDatabaseEditDate"), "sha256": db.sha256,
            "source": meta.get("url"), "createdBy": meta.get("createdBy"),
            "maintainedBy": meta.get("maintainedBy"), "copyright": meta.get("copyright"),
            "pdbChecked": online,
            # Distinct structures. Every structure that is not a parent counts as a
            # daughter structure: mpstruc's <memberProtein> records and its related
            # PDB entries alike. Cross-listed structures are counted once.
            "counts": _counts(db),
            "flagCounts": {s: sum(1 for f in findings if f["severity"] == s) for s in ("error", "warning", "info")},
        },
        "groups": groups,
        "relPdb": rel_info,
    }


# ---- EMDB links, patched into the published files in place -------------------
_JSON_STYLES = [  # how this project writes its data files; the first that reproduces a file is kept
    {"ensure_ascii": False, "separators": (",", ":")},
    {"separators": (",", ":")},
    {"indent": 1, "ensure_ascii": False},
    {"indent": 1},
    {"indent": 2, "ensure_ascii": False},
    {"indent": 2},
]


def _json_style(text: str, doc) -> dict | None:
    body = text[:-1] if text.endswith("\n") else text
    for st in _JSON_STYLES:
        if json.dumps(doc, **st) == body:
            return st
    return None


def _set_emdb(rec: dict, ids: list[str] | None) -> bool:
    """Set rec["emdb"] to ids (dropped when there are none); True if rec changed."""
    old = rec.get("emdb")
    if ids:
        if old == ids:
            return False
        rec["emdb"] = list(ids)
        return True
    if "emdb" in rec:
        del rec["emdb"]
        return True
    return False


def _mp_records(doc: dict):
    """(PDB ID, record to hold "emdb") for every parent, daughter and related structure."""
    for g in doc.get("groups") or []:
        for sg in g.get("subgroups") or []:
            for p in sg.get("parents") or []:
                for e in [p, *(p.get("d") or [])]:
                    if isinstance(e.get("x"), dict) and e.get("pdb"):
                        yield str(e["pdb"]).upper(), e["x"]
    for pid, info in (doc.get("relPdb") or {}).items():
        if isinstance(info, dict):
            yield str(pid).upper(), info


def _cand_records(doc: dict):
    """(PDB ID, candidate) for every suggestion, same-paper siblings included."""
    for c in doc.get("candidates") or []:
        if isinstance(c, dict) and c.get("pdb"):
            yield str(c["pdb"]).upper(), c


def patch_links(mp_path: Path = OUT, cand_path: Path | None = CANDIDATES, fetch=None) -> dict:
    """Refresh "emdb" (EMDB map IDs) on every structure of the existing site files.

    data/mpstruc.json: each parent's and daughter's "x" and each "relPdb" record;
    data/candidates.json: each candidate (siblings are candidates too). Nothing
    else is touched; each file keeps its own JSON formatting and is only rewritten
    when something changed, so running it twice changes nothing the second time.
    fetch: rcsb.emdb_ids, or a stand-in for tests. Returns counts per file.
    """
    if fetch is None:
        from .rcsb import emdb_ids as fetch
    files = [(Path(mp_path), _mp_records)]
    if cand_path is not None and Path(cand_path).exists():
        files.append((Path(cand_path), _cand_records))
    loaded = []
    for path, records in files:
        text = path.read_text(encoding="utf-8")
        doc = json.loads(text)
        style = _json_style(text, doc)
        if style is None:
            raise SystemExit(f"{path}: unknown JSON formatting; not rewritten")
        loaded.append((path, records, text, doc, style))
    ids = sorted({pid for _, records, _, doc, _ in loaded for pid, _ in records(doc)})
    maps = {str(k).upper(): v for k, v in fetch(ids).items()}
    out = {}
    for path, records, text, doc, style in loaded:
        changed, with_maps, n = 0, set(), set()
        for pid, rec in records(doc):
            got = [m for m in maps.get(pid) or [] if re.fullmatch(r"EMD-\d{4,6}", str(m))]
            changed += _set_emdb(rec, got)
            n.add(pid)
            if got:
                with_maps.add(pid)
        if changed:
            new = json.dumps(doc, **style) + ("\n" if text.endswith("\n") else "")
            path.write_text(new, encoding="utf-8")
        out[path.name] = {"structures": len(n), "withEmdb": len(with_maps), "recordsChanged": changed}
    return out


def main(argv: list[str] | None = None) -> None:
    import sys
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["links"]:
        lp = argparse.ArgumentParser(prog="python -m mpstruc.web links",
                                     description="Refresh the EMDB map IDs in the existing site data files, in place.")
        lp.add_argument("--mpstruc", default=str(OUT))
        lp.add_argument("--candidates", default=str(CANDIDATES))
        la = lp.parse_args(argv[1:])
        for name, c in patch_links(Path(la.mpstruc), Path(la.candidates)).items():
            print(f"{name}: {c['withEmdb']} of {c['structures']} structures have EMDB maps; "
                  f"{c['recordsChanged']} records updated")
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)
    data = build(load(source_xml()), online=not args.offline)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")))
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
