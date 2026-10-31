"""Data-integrity checks over a parsed mpstruc Database.

Each check returns a list of findings. Nothing here changes the data: the
report is for curators to review. Some findings are intentional curation
(cross-listings) and are labelled as such rather than as errors.
"""
from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field

from .model import Database, Entry
from .resolution import parse_resolution
from .text import (ENTITY_REPAIRS, broken_entities, mentioned_pdb_ids, to_text,
                   unknown_entities, unterminated_entities)

PDB_ID_FORMAT = re.compile(r"[0-9][A-Z0-9]{3}\Z")
KNOWN_DOMAINS = {"Eukaryota", "Bacteria", "Archaea", "Viruses"}
ALLOWED_TAGS = {"br", "sub", "sup", "em", "i", "b", "a"}
TAG_NAME = re.compile(r"<\s*/?\s*([A-Za-z][A-Za-z0-9]*)")


@dataclass
class Finding:
    check: str
    severity: str  # "error" | "warning" | "info"
    key: str | None
    pdb: str | None
    message: str
    data: dict = field(default_factory=dict)


def _fields(e: Entry):
    yield "name", e.name
    yield "species", e.species
    yield "expressedInSpecies", e.expressedInSpecies
    yield "resolution", e.resolution
    yield "description", e.description
    for i, b in enumerate([e.bibliography, *e.secondary]):
        for f in ("authors", "title", "journal", "notes"):
            yield f"bibliography[{i}].{f}", getattr(b, f)


def check_structure(db: Database) -> list[Finding]:
    out = []
    names = Counter(sg.name for sg in db.subgroups())
    for sg in db.subgroups():
        if not sg.parents:
            out.append(Finding("empty_subgroup", "info", sg.key, None,
                               f"Subgroup '{to_text(sg.name)}' has no entries (header subgroup?)"))
        if names[sg.name] > 1:
            out.append(Finding("duplicate_subgroup_name", "info", sg.key, None,
                               f"Subgroup name '{to_text(sg.name)}' is used by {names[sg.name]} subgroups; "
                               "identify subgroups by key"))
        for p in sg.parents:
            for d in p.daughters:
                if d.masterProteinPdbCode != p.pdbCode:
                    out.append(Finding("master_code_mismatch", "error", d.key, d.pdbCode,
                                       f"masterProteinPdbCode {d.masterProteinPdbCode} != enclosing parent {p.pdbCode}"))
    for e in db.entries():
        ids = [e.pdbCode, *[r.pdb for r in e.related]]
        for pid in ids:
            if not PDB_ID_FORMAT.match(pid.upper()):
                out.append(Finding("bad_pdb_id", "error", e.key, pid, f"'{pid}' is not a valid PDB ID"))
        dup = [k for k, n in Counter(r.pdb for r in e.related).items() if n > 1]
        if dup:
            out.append(Finding("duplicate_related", "warning", e.key, e.pdbCode,
                               f"Related list repeats {', '.join(dup)}"))
        if e.pdbCode in {r.pdb for r in e.related}:
            out.append(Finding("self_related", "warning", e.key, e.pdbCode, "Entry lists itself as related"))
    return out


def check_cross_listings(db: Database) -> list[Finding]:
    where: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    sg_of = {}
    for g in db.groups:
        for sg in g.subgroups:
            for p in sg.parents:
                for e in (p, *p.daughters):
                    sg_of[e.key] = sg
                    where[e.pdbCode].append((e.key, e.kind, sg.key))
                    for r in e.related:
                        where[r.pdb].append((r.key, "related", sg.key))
    out = []
    for pid, occ in where.items():
        if len(occ) < 2:
            continue
        sgs = {s for _, _, s in occ}
        documented = any("listed under" in (to_text(e.description) or "").lower()
                         for e in db.entries() if e.pdbCode == pid)
        if len(sgs) > 1:
            out.append(Finding("cross_listing", "info", occ[0][0], pid,
                               f"Listed in {len(sgs)} subgroups" + (" (noted in description)" if documented else ""),
                               {"occurrences": occ, "documented": documented}))
        else:
            # An entry listing itself (reported as self_related) or its own
            # parent/daughter as related is a family-internal link, not a repeat.
            roots = {k.split("-D")[0].split("-R")[0] for k, _, _ in occ}
            if len(roots) == 1:
                continue
            out.append(Finding("repeated_in_subgroup", "warning", occ[0][0], pid,
                               f"Listed {len(occ)} times within subgroup {occ[0][2]}: "
                               + ", ".join(k for k, _, _ in occ), {"occurrences": occ}))
    return out


def check_markup(db: Database) -> list[Finding]:
    out = []
    for e in db.entries():
        for fname, raw in _fields(e):
            if not raw:
                continue
            for name in unknown_entities(raw):
                fix = ENTITY_REPAIRS.get(name)
                out.append(Finding("unknown_entity", "warning", e.key, e.pdbCode,
                                   f"{fname}: &{name}; is not an HTML entity"
                                   + (f" (displayed as '{fix}')" if fix else " (left as text)")))
            for ent in broken_entities(raw):
                out.append(Finding("broken_entity", "warning", e.key, e.pdbCode,
                                   f"{fname}: malformed '{ent}'"
                                   + (" (displayed as 'Å')" if "Aring" in ent else " (left as '¨'; probably ü)")))
            if re.search("[\u200b\u200c\u200d\u2060\ufeff]", raw):
                out.append(Finding("invisible_characters", "info", e.key, e.pdbCode,
                                   f"{fname}: contains zero-width characters (removed for display)"))
            for ent in unterminated_entities(raw):
                out.append(Finding("unterminated_entity", "warning", e.key, e.pdbCode,
                                   f"{fname}: '{ent}' is missing its ';'"))
            if re.search(r"&(?![A-Za-z][A-Za-z0-9]*;?|#)", raw):
                out.append(Finding("bare_ampersand", "info", e.key, e.pdbCode,
                                   f"{fname}: unescaped '&' (invalid XML)"))
            tags = [t.lower() for t in TAG_NAME.findall(raw)]
            bad = sorted({t for t in tags if t not in ALLOWED_TAGS})
            if bad:
                out.append(Finding("malformed_markup", "warning", e.key, e.pdbCode,
                                   f"{fname}: unexpected tag(s) {', '.join('<' + t + '>' for t in bad)}"))
            opens = len(re.findall(r"<a[\s>]", raw, re.I))
            closes = len(re.findall(r"</a\s*>", raw, re.I))
            if opens != closes:
                out.append(Finding("unclosed_link", "warning", e.key, e.pdbCode,
                                   f"{fname}: {opens} <a> vs {closes} </a>"))
            if re.search(r"<(sub|sup|em|i|b)>\s*</\1>", raw):
                out.append(Finding("empty_tag", "info", e.key, e.pdbCode, f"{fname}: empty formatting tag"))
    return out


def check_fields(db: Database) -> list[Finding]:
    out = []
    for e in db.entries():
        for f in ("name", "species", "taxonomicDomain"):
            if not to_text(getattr(e, f)):
                out.append(Finding("missing_field", "warning", e.key, e.pdbCode, f"{f} is empty"))
        dom = to_text(e.taxonomicDomain)
        if dom and dom not in KNOWN_DOMAINS:
            out.append(Finding("unknown_domain", "warning", e.key, e.pdbCode, f"taxonomicDomain '{dom}'"))
        res = parse_resolution(e.resolution)
        if res.kind == "unparseable":
            out.append(Finding("resolution_unparseable", "warning", e.key, e.pdbCode,
                               f"resolution '{to_text(e.resolution)}' is neither a number nor a method"))
        elif res.kind == "missing":
            out.append(Finding("resolution_missing", "warning", e.key, e.pdbCode, "resolution is empty"))
        elif res.kind == "multiple":
            out.append(Finding("resolution_multiple", "info", e.key, e.pdbCode,
                               f"resolution '{to_text(e.resolution)}' has several values; using {res.value}"))
        for i, b in enumerate([e.bibliography, *e.secondary]):
            label = "primary" if i == 0 else f"secondary {i}"
            if not to_text(b.title):
                out.append(Finding("citation_incomplete", "warning", e.key, e.pdbCode, f"{label} citation has no title"))
            if not to_text(b.pubMedId) and not to_text(b.doi):
                out.append(Finding("citation_unlinked", "info", e.key, e.pdbCode,
                                   f"{label} citation has neither PubMed ID nor DOI"))
            pm = to_text(b.pubMedId)
            if pm and not pm.isdigit():
                out.append(Finding("bad_pubmed_id", "warning", e.key, e.pdbCode, f"{label} PubMed ID '{pm}'"))
            doi = to_text(b.doi)
            if doi and not re.fullmatch(r"10\.\d{4,9}/\S+", doi):
                out.append(Finding("bad_doi", "warning", e.key, e.pdbCode, f"{label} DOI '{doi}'"))
            yr = to_text(b.year)
            if yr and not (yr.isdigit() and 1950 <= int(yr) <= 2100):
                out.append(Finding("bad_year", "warning", e.key, e.pdbCode, f"{label} year '{yr}'"))
    return out


PREPRINT_DOI_PREFIXES = ("10.1101/", "10.21203/", "10.64898/", "10.48550/", "10.20944/")
PREPRINT_SOURCES = {"biorxiv", "res sq", "medrxiv", "arxiv", "preprints"}
GENERIC_WORDS = {
    "structure", "structures", "protein", "proteins", "complex", "crystal", "with", "bound", "from",
    "state", "states", "human", "mutant", "form", "cryo", "domain", "in", "the", "and", "of", "bacterial",
    "membrane", "transmembrane", "basis", "structural", "insights", "into", "mechanism", "reveals",
    "resolution", "high", "analysis", "function", "molecular", "determination", "complexes", "family",
}


def _norm_doi(d) -> str:
    d = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:\s*)", "", str(d or "").strip().lower())
    return d.rstrip(".")


def _pmid(v) -> str:
    """PubMed ID as a string, or '' for missing values and RCSB's -1 sentinel."""
    s = str(v or "").strip()
    return s if s.isdigit() and int(s) > 0 else ""


def _sim(a: str | None, b: str | None) -> float:
    from difflib import SequenceMatcher
    na = re.sub(r"[^a-z0-9]+", " ", (a or "").lower()).strip()
    nb = re.sub(r"[^a-z0-9]+", " ", (b or "").lower()).strip()
    return SequenceMatcher(None, na, nb).ratio() if na and nb else 0.0


GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon", "κ": "kappa",
         "μ": "mu", "σ": "sigma", "ζ": "zeta", "φ": "phi", "ω": "omega", "θ": "theta"}


def _content_words(s: str | None) -> set[str]:
    s = (s or "").lower()
    for g, name in GREEK.items():
        s = s.replace(g, name)
    words = set(re.findall(r"[a-z0-9]+", s))
    return {w for w in words if w not in GENERIC_WORDS and not w.isdigit()
            and (len(w) >= 3 or any(c.isdigit() for c in w))}


def _ids(doi: str, pmid: str, pm: dict) -> tuple[set[str], set[str]]:
    """All DOIs and PMIDs a citation resolves to (PubMed adds the DOI of a PMID)."""
    dois = {d for d in (doi, _norm_doi((pm.get(pmid) or {}).get("doi"))) if d}
    return dois, ({pmid} if pmid else set())


def _same_paper(a, b) -> bool:
    (ad, ap), (bd, bp) = a, b
    squash = lambda ds: {d.replace(".", "") for d in ds}  # tolerate a missing '.' (malformed DOIs)
    return bool(ad & bd or ap & bp or squash(ad) & squash(bd))


def _is_preprint(doi: str, pmid: str, pm: dict) -> bool:
    src = ((pm.get(pmid) or {}).get("source") or "").lower()
    return doi.startswith(PREPRINT_DOI_PREFIXES) or src in PREPRINT_SOURCES


def check_against_pdb(db: Database) -> list[Finding]:
    """Compare curated citations, resolutions and identity with the PDB and PubMed.

    Needs network (RCSB GraphQL, NCBI esummary); cached under cache/. Two
    citations are the same paper if any DOI or PubMed ID agrees after
    resolving PMIDs to DOIs, so a DOI on one side and a PMID on the other
    is not reported. Nothing is reported when the PDB has no citation IDs.
    """
    from . import pubmed, rcsb
    ents = list(db.entries())
    pdb = rcsb.entry_summaries([e.pdbCode for e in ents])
    pmids = [_pmid(to_text(b.pubMedId)) for e in ents for b in (e.bibliography, *e.secondary)]
    pmids += [_pmid(((s or {}).get("rcsb_primary_citation") or {}).get("pdbx_database_id_PubMed"))
              for s in pdb.values()]
    pm = pubmed.summaries([p for p in pmids if p])
    out = []
    for e in ents:
        # PubMed ID vs the curated title, for every bibliography record.
        for i, b in enumerate([e.bibliography, *e.secondary]):
            pid = _pmid(to_text(b.pubMedId))
            rec = pm.get(pid) if pid else None
            title = to_text(b.title)
            if pid and rec is None:
                out.append(Finding("pubmed_id_not_found", "warning", e.key, e.pdbCode,
                                   f"{'primary' if i == 0 else 'secondary'} PubMed ID {pid} does not exist in PubMed"))
            elif rec and title and _sim(rec["title"], title) < 0.6:
                out.append(Finding("pubmed_title_mismatch", "warning", e.key, e.pdbCode,
                                   f"{'Primary' if i == 0 else 'Secondary'} PubMed ID {pid} is "
                                   f"'{rec['title']}' ({rec['source']} {rec['year']}), but the citation title is '{title}'",
                                   {"pubmed": pid, "pubmed_title": rec["title"]}))

        s = pdb.get(e.pdbCode)
        if not s:
            continue
        cit = s.get("rcsb_primary_citation") or {}
        p_doi, p_pm = _norm_doi(cit.get("pdbx_database_id_DOI")), _pmid(cit.get("pdbx_database_id_PubMed"))
        bibs = [e.bibliography, *e.secondary]
        ours = [(_norm_doi(to_text(b.doi)), _pmid(to_text(b.pubMedId))) for b in bibs]
        theirs = _ids(p_doi, p_pm, pm)
        o_doi, o_pm = ours[0]
        pdb_title = (s.get("struct") or {}).get("title")
        methods = [m.get("method", "") for m in s.get("exptl") or []]

        # A PDB ID that names an unrelated structure: nothing in the curated
        # name appears in the PDB's title or citation title.
        name_words = _content_words(to_text(e.name))
        pdb_words = _content_words(pdb_title) | _content_words(cit.get("title"))
        unrelated = bool(name_words) and not (name_words & pdb_words)

        if (p_doi or p_pm) and (o_doi or o_pm):
            primary = _ids(o_doi, o_pm, pm)
            if _same_paper(primary, theirs):
                if o_pm and p_pm and o_pm != p_pm and not (_ids("", o_pm, pm)[0] & theirs[0]):
                    out.append(Finding("pubmed_id_conflict", "warning", e.key, e.pdbCode,
                                       f"DOI matches the PDB, but PubMed ID {o_pm} differs from the PDB's {p_pm}"
                                       + (f" ({pm[o_pm]['title']})" if pm.get(o_pm) else "")))
                elif o_doi and p_doi and o_doi != p_doi and o_doi.replace(".", "") != p_doi.replace(".", ""):
                    pm_doi = _norm_doi((pm.get(o_pm or p_pm) or {}).get("doi"))
                    if pm_doi and pm_doi == o_doi:
                        out.append(Finding("doi_conflict", "info", e.key, e.pdbCode,
                                           f"The PDB's DOI '{p_doi}' looks wrong; mpstruc's '{o_doi}' agrees with PubMed"))
                    else:
                        out.append(Finding("doi_conflict", "warning", e.key, e.pdbCode,
                                           f"PubMed ID matches the PDB, but DOI '{o_doi}' differs from the PDB's '{p_doi}'"
                                           + (" (truncated?)" if p_doi.startswith(o_doi) else "")
                                           + (f"; PubMed has '{pm_doi}'" if pm_doi else "")))
            elif any(_same_paper(_ids(d, p, pm), theirs) for d, p in ours[1:]):
                out.append(Finding("citation_differs_from_pdb", "info", e.key, e.pdbCode,
                                   "PDB primary citation is listed here as a secondary reference"))
            elif unrelated:
                out.append(Finding("pdb_id_points_elsewhere", "error", e.key, e.pdbCode,
                                   f"PDB {e.pdbCode} is '{pdb_title}' ({', '.join(methods)}), which does not match "
                                   f"'{to_text(e.name)}'; its citation differs too. Wrong PDB ID?",
                                   {"pdb_title": pdb_title}))
            else:
                ours_title, their_title = to_text(e.bibliography.title), cit.get("title")
                # Preprints are often retitled on publication, so a looser
                # threshold applies when one side is a preprint.
                sim = _sim(ours_title, their_title)
                ours_pre, theirs_pre = _is_preprint(o_doi, o_pm, pm), _is_preprint(p_doi, p_pm, pm)
                same_work = sim >= 0.6 if (ours_pre or theirs_pre) else sim >= 0.95
                if same_work and ours_pre and not theirs_pre:
                    out.append(Finding("citation_preprint", "warning", e.key, e.pdbCode,
                                       f"mpstruc cites the preprint; the published version is "
                                       f"{cit.get('journal_abbrev')} {cit.get('year')} (DOI {p_doi or '—'}, PMID {p_pm or '—'})"))
                elif same_work and theirs_pre:
                    pass  # mpstruc already cites the published paper; the PDB still has the preprint
                elif same_work:
                    out.append(Finding("citation_ids_wrong", "warning", e.key, e.pdbCode,
                                       f"Same title as the PDB's citation, but the PubMed ID/DOI differ "
                                       f"(mpstruc {o_pm or '—'} / {o_doi or '—'}, PDB {p_pm or '—'} / {p_doi or '—'})"))
                else:
                    out.append(Finding(
                        "citation_differs_from_pdb", "warning", e.key, e.pdbCode,
                        f"Primary citation '{ours_title}' differs from the PDB's '{their_title}' ({cit.get('year')})",
                        {"pdb_doi": p_doi, "pdb_pubmed": p_pm, "mpstruc_doi": o_doi, "mpstruc_pubmed": o_pm}))
        elif unrelated and pdb_title:
            out.append(Finding("pdb_title_mismatch", "warning", e.key, e.pdbCode,
                               f"PDB {e.pdbCode} is titled '{pdb_title}', which shares no words with '{to_text(e.name)}'"))

        res = parse_resolution(e.resolution)
        pres = (s.get("rcsb_entry_info") or {}).get("resolution_combined") or []
        nmr_only = bool(methods) and all("NMR" in m for m in methods)
        if res.value is not None and nmr_only:
            out.append(Finding("resolution_vs_method", "warning", e.key, e.pdbCode,
                               f"mpstruc gives {res.value} Å, but the PDB entry is {', '.join(methods)} (no resolution)"))
        elif res.method_hint and "NMR" in res.method_hint and methods and not any("NMR" in m for m in methods):
            out.append(Finding("resolution_vs_method", "warning", e.key, e.pdbCode,
                               f"mpstruc says '{to_text(e.resolution)}', but the PDB entry is {', '.join(methods)}"
                               + (f" at {', '.join(map(str, pres))} Å" if pres else "")))
        elif res.value is not None and pres and min(abs(res.value - r) for r in pres) > 0.051:
            out.append(Finding("resolution_differs_from_pdb", "warning", e.key, e.pdbCode,
                               f"mpstruc resolution {res.value} Å, PDB reports {', '.join(map(str, pres))} Å"))
    return out


def check_mentions(db: Database) -> list[Finding]:
    """mpstruc descriptions cite each related structure ('..., 2.22 Å: 6BL3').
    Compare those citations with the entry's own related list, both ways."""
    where: dict[str, list[str]] = defaultdict(list)
    for e in db.entries():
        where[e.pdbCode.upper()].append(e.key)
        for r in e.related:
            where[r.pdb.upper()].append(r.key)
    out = []
    for e in db.entries():
        mentioned = mentioned_pdb_ids(e.description)
        own = {e.pdbCode.upper()} | {r.pdb.upper() for r in e.related}
        missing = [m for m in mentioned if m not in own]
        if missing:
            elsewhere = {m: where[m] for m in missing if where.get(m)}
            out.append(Finding("described_not_related", "warning", e.key, e.pdbCode,
                               f"Description cites {', '.join(missing)} but its related list does not include "
                               + ("them" if len(missing) > 1 else "it")
                               + (" (listed elsewhere: " + "; ".join(f"{m} at {', '.join(k)}" for m, k in elsewhere.items()) + ")"
                                  if elsewhere else ""),
                               {"pdb_ids": missing, "elsewhere": elsewhere}))
        undescribed = [r.pdb for r in e.related if r.pdb.upper() not in mentioned]
        if undescribed and mentioned:
            out.append(Finding("related_not_described", "info", e.key, e.pdbCode,
                               f"Related {', '.join(undescribed)} not cited in the description",
                               {"pdb_ids": undescribed}))
    return out


def check_obsolete(db: Database) -> list[Finding]:
    """Needs network (RCSB holdings); results are cached under cache/."""
    from . import rcsb
    current = rcsb.current_entry_ids()
    unreleased = rcsb.unreleased_entry_ids()
    occ: dict[str, str] = {}
    for e in db.entries():
        occ.setdefault(e.pdbCode, e.key)
        for r in e.related:
            occ.setdefault(r.pdb, r.key)
    out = []
    for pid, key in occ.items():
        if pid in current:
            continue
        if pid.upper() in current:
            out.append(Finding("pdb_id_case", "warning", key, pid,
                               f"'{pid}' should be '{pid.upper()}' (PDB IDs are upper case)"))
            continue
        if pid in unreleased:
            out.append(Finding("pdb_unreleased", "info", key, pid, "Not yet released by the PDB"))
            continue
        rem = rcsb.removed_entry(pid)
        if rem:
            by = rem.get("id_codes_replaced_by") or []
            msg = f"superseded: {pid} → {', '.join(by)}" if by else f"obsolete: {pid} (no replacement)"
            out.append(Finding("pdb_obsolete", "warning", key, pid, msg,
                               {"replaced_by": by, "removed": rem.get("remove_date")}))
        else:
            out.append(Finding("pdb_unknown", "error", key, pid, "Not a current, unreleased or removed PDB entry"))
    return out


OFFLINE_CHECKS = [check_structure, check_cross_listings, check_markup, check_fields, check_mentions]
ONLINE_CHECKS = [check_obsolete, check_against_pdb]


def run(db: Database, online: bool = False) -> list[dict]:
    findings = [f for chk in OFFLINE_CHECKS for f in chk(db)]
    if online:
        findings += [f for chk in ONLINE_CHECKS for f in chk(db)]
    return [asdict(f) for f in findings]
