"""Draft mpstruc entries for new PDB structures.

Builds the fields of an mpstruc <protein>/<memberProtein> record (name,
species, taxonomic domain, expression host, resolution, description and
bibliography) from the PDB, written the way mpstruc writes them:
the protein named as earlier mpstruc entries for the same protein name it,
partners and ligands in the name ("Dopamine D<sub>1</sub> receptor -
G<sub>s</sub> complex with bound SKF-83566"), "cryo-EM structure" as the
description, related structures as "with bound X, 3.10 Å: 9ABC" (only what
differs from the entry they are listed under), and authors as
"Ahmed S, Huang M, & Lee CH".
"""
from __future__ import annotations

import html
import re
from collections import Counter

DOMAINS = ("Eukaryota", "Bacteria", "Archaea", "Viruses")
# Ions written the way mpstruc writes them.
IONS = {"CA": "Ca<sup>2+</sup>", "MG": "Mg<sup>2+</sup>", "ZN": "Zn<sup>2+</sup>", "MN": "Mn<sup>2+</sup>",
        "NA": "Na<sup>+</sup>", "K": "K<sup>+</sup>", "CL": "Cl<sup>-</sup>", "CU": "Cu<sup>2+</sup>",
        "CU1": "Cu<sup>+</sup>", "FE": "Fe<sup>3+</sup>", "FE2": "Fe<sup>2+</sup>", "CD": "Cd<sup>2+</sup>",
        "LI": "Li<sup>+</sup>", "RB": "Rb<sup>+</sup>", "CS": "Cs<sup>+</sup>", "BA": "Ba<sup>2+</sup>",
        "SR": "Sr<sup>2+</sup>", "NI": "Ni<sup>2+</sup>", "CO": "Co<sup>2+</sup>", "TL": "Tl<sup>+</sup>"}
# Water, buffers, cryoprotectants, detergents and common lipids: never "bound ligands" in mpstruc.
ADDITIVES = set("""HOH DOD GOL EDO PEG PG4 PGE 1PE P6G PE4 PE5 SO4 PO4 ACT FMT DMS MPD TRS EPE MES BME IMD
CIT TLA ACY IPA EOH MOH NO3 IOD BR UNX UNL BOG LMT LMN LDA C8E BNG HTG OLC OLA OLB 78M 78N 7E9 MPG
Y01 CLR PCW POV PEE PLC LPP PGW CDL 3PE PC1 6PL P5S PEF PEV D10 D12 OCT HEX DD9 R16 TRD LFA PLM MYR
STE PA8 PA6 PTY NAG NDG MAN BMA FUC GAL GLC SIA DMU UMQ CPS SDS DPC LDM 9Z9 AV0 PX4 PCF PGT U10
7PH PTD 8PE 44E 2CV 17F J40 P0E LHG LMG DGD SQD""".split())
ANTIBODY = re.compile(r"\b(fab|nanobody|nb\d*|scfv|antibody|heavy chain|light chain|sybody|darpin|megabody|"
                      r"igg|vhh|monobody)\b", re.I)


def _authors(names: list[str] | None) -> str | None:
    """['Ahmed, S.', 'Lee, C.H.'] -> 'Ahmed S, & Lee CH' (mpstruc style)."""
    out = []
    for n in names or []:
        last, _, first = n.partition(",")
        initials = re.sub(r"[^A-Za-z]", "", first)
        out.append(f"{last.strip()} {initials}".strip())
    if not out:
        return None
    return out[0] if len(out) == 1 else ", ".join(out[:-1]) + ", & " + out[-1]


def _cap(s: str) -> str:
    return s[:1].upper() + s[1:] if s else s


FUSION_NAME = re.compile(r"^(green fluorescent protein|gfp|egfp|t4 lysozyme|endolysin|lysozyme|soluble cytochrome b562|"
                         r"bril|maltose[- ]binding.*|maltodextrin[- ]binding.*|thioredoxin.*|rubredoxin|"
                         r"glycogen synthase.*|mcherry|venus)$", re.I)


def _clean_name(desc: str) -> str:
    """PDB entity name -> mpstruc-like protein name: drop fusion partners
    ('X,Green fluorescent protein') and 'Isoform 1 of'."""
    parts = [p.strip() for p in desc.split(",") if p.strip()]
    kept = [p for p in parts if not FUSION_NAME.match(p)] or parts
    name = ", ".join(kept)
    return re.sub(r"^isoform \w+ of ", "", name, flags=re.I)


def _method(methods: list[str]) -> str | None:
    m = " ".join(methods).upper()
    if "ELECTRON MICROSCOPY" in m:
        return "cryo-EM structure"
    if "NMR" in m:
        return "solid-state NMR structure" if "SOLID" in m else "NMR structure"
    if "ELECTRON CRYSTALLOGRAPHY" in m:
        return "electron crystallography structure"
    if "X-RAY" in m:
        return "x-ray structure"
    return None


# Words in a PDB title that are never the name of a bound ligand.
NOT_LIGAND = set("""the a an and or of with without in to by from for its apo ligand ligands agonist agonists antagonist
inhibitor inhibitors drug compound complex complexes state form structure receptor channel transporter protein
signaling signalling active inactive bound unbound full partial inverse allosteric modulator positive negative
nanobody fab scfv antibody nanodisc nanodiscs lipid detergent cryo-em crystal""".split())
TITLE_LIGAND = re.compile(r"\b(?!with\b)([A-Za-z0-9][\w()',+\-]*?)(?=[- ]bound\b(?! to))"
                          r"|\b(?:bound to|with(?: bound)?|by|and)\s+(?:the\s+)?(?:(?:full|partial|inverse|positive|negative)\s+)?"
                          r"(?:(?:agonist|antagonist|inhibitor|ligand|drug|allosteric modulator)\s+)?"
                          r"((?:compound|cmpd)\s+[\w\-]+|[A-Za-z0-9][\w()',+\-]{2,})")
# Stabilisers mpstruc leaves out of G protein complexes.
G_STABILISER = re.compile(r"^(nb ?35|nanobody[- ]?35|scfv ?16|single[- ]chain.*16|fab ?16)$", re.I)
PEPTIDE_MAX = 40  # partner chains this short are peptide ligands ("with bound C5a")


def _title_ligands(title: str | None, entity_names: list[str]) -> list[str]:
    """Ligand names the PDB title gives: 'SKF-83566 bound dopamine receptor' -> ['SKF-83566']."""
    words = {w.lower() for n in entity_names for w in re.split(r"[\s,/]+", n) if w}
    out = []
    for m in TITLE_LIGAND.finditer(title or ""):
        tok = (m.group(1) or m.group(2) or "").strip("-,'()")
        low = tok.lower()
        if (len(tok) < 3 or low in NOT_LIGAND or low in words or re.fullmatch(r"\d+(\.\d+)?", tok)
                or re.search(r"(?i)^(g[sioqt]|minig\w*|nb\d+|scfv\d*|fab\w*)$|-g[sioq]$", tok)
                or any(low in w or w in low for w in words if len(w) > 3 and re.search(r"\d", w))):
            continue
        if re.fullmatch(r"[A-Z][a-z]+", tok):
            tok = low  # 'Histamine-bound ...' -> histamine
        if tok not in out:
            out.append(tok)
    return out


def _ligands(det: dict, entity_names: list[str] = ()) -> list[str]:
    """Bound ligands as mpstruc names them: ions as Ca<sup>2+</sup>, drugs by
    their common name. A systematic chemical name is replaced by the name the
    PDB title uses ('SKF-83566'), and only kept as the 3-character code when
    the title gives none."""
    out, coded = [], []
    for n in det.get("nonpolymer_entities") or []:
        comp = ((n.get("nonpolymer_comp") or {}).get("chem_comp") or {})
        cid, name = comp.get("id") or "", comp.get("name") or ""
        if not cid or cid in ADDITIVES:
            continue
        if cid in IONS:
            label = IONS[cid]
        elif len(name) <= 24 and not re.search(r"\d,\d|\(\d|~|\{", name):
            label = name.lower() if name.isupper() else name
        else:
            label = cid
            coded.append(cid)
        if label not in out:
            out.append(label)
    if coded:
        named = [t for t in _title_ligands((det.get("struct") or {}).get("title"), entity_names)
                 if t.lower() not in {o.lower() for o in out}]
        for cid, t in zip(coded, named):
            out[out.index(cid)] = t
    return out[:4]


def _join(items: list[str]) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _gprotein(alpha: list[str], mini: bool) -> str:
    """'Guanine nucleotide-binding protein G(s) subunit alpha ...' -> 'G<sub>s</sub>'."""
    m = next((re.search(r"G\((\w+)\)(?:,? subunit alpha| alpha)", a, re.I) for a in alpha), None)
    sub = m.group(1).lower() if m else ""
    sub = {"i": "i", "s": "s", "o": "o", "q": "q", "t": "t", "olf": "olf", "z": "z"}.get(sub, sub.upper() if sub.isdigit() else sub)
    if not sub:
        return "G protein"
    return f"{'miniG' if mini else 'G'}<sub>{sub}</sub>"


def parts(det: dict, core: set[str]) -> dict:
    """What a structure contains besides its membrane protein(s): the G protein
    (if any), other partner proteins, and bound ligands (peptides included)."""
    partners, peptides, alpha, names, gbg = [], [], [], [], False
    for e in det.get("polymer_entities") or []:
        poly = e.get("entity_poly") or {}
        if poly.get("rcsb_entity_polymer_type") != "Protein":
            continue
        ups = set((e.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or [])
        desc = ((e.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()
        names.append(desc)
        if not desc or ups & core:
            continue
        if re.match(r"(engineered |mini[- ]?)?guanine nucleotide-binding protein|mini[- ]?g", desc, re.I):
            if re.search(r"alpha|G\(\w+\) subunit$|mini", desc, re.I):
                alpha.append(desc)
            gbg = True
            continue
        n = poly.get("rcsb_sample_sequence_length")
        target = peptides if n and n <= PEPTIDE_MAX and not ANTIBODY.search(desc) else partners
        if desc not in target:
            target.append(desc)
    title = (det.get("struct") or {}).get("title") or ""
    g = None
    if gbg:
        g = _gprotein(alpha, bool(re.search(r"mini[- ]?g", title + " ".join(alpha), re.I)))
        partners = [p for p in partners if not G_STABILISER.match(p)]
    return {"g": g, "partners": partners[:3], "ligands": peptides + _ligands(det, names)}


def phrase(pt: dict, base: dict | None = None) -> str:
    """mpstruc wording for a structure's contents: '- G<sub>s</sub> complex with
    bound X' / 'in complex with Y, with bound X'. With `base`, only what
    differs from it (for a related-structure line)."""
    g, partners, ligs = pt["g"], pt["partners"], pt["ligands"]
    if base is not None:
        g = g if g != base["g"] else None
        partners = [p for p in partners if p not in base["partners"]]
        if ligs == base["ligands"]:
            ligs = []
        elif not ligs and base["ligands"]:
            return ", ".join(x for x in (f"{g} complex" if g else "", "apo form") if x)
    out = f"{'- ' if base is None else ''}{g} complex" if g else ""
    if partners:
        out = f"{out} in complex with {_join(partners)}".strip()
    if ligs:
        out = f"{out}{',' if partners else ''} with bound {_join(ligs)}".strip()
    return out


def condition(det: dict, core: set[str]) -> str:
    """The phrase for one structure on its own, e.g. 'in complex with X, with bound Y'."""
    return phrase(parts(det, core)).removeprefix("- ")


CUT = re.compile(r"\s*-\s*(?:mini\s*)?G(?:<sub>|\s*protein\b|\s*[sioq]\b)"
                 r"|\s+(?:complexed|in complex|bound|together)\s+(?:with|to)\b|\s+with\s|,\s"
                 r"|\s+in\s+(?=(?:lipid|nanodisc|detergent|saposin|amphipol|peptidisc|liposome|the|a|an|complex"
                 r"|[A-Z0-9]{3,})\b)|\s+(?:apo|bound|unbound|ligand-free|inward|outward|open|closed)\b", re.I)


def stem(name: str) -> str:
    """An mpstruc name without its condition: 'Dopamine D<sub>1</sub> receptor
    complexed with G<sub>s</sub> protein with bound SKF-81297' -> 'Dopamine
    D<sub>1</sub> receptor'. Text in parentheses is never cut."""
    depth = 0
    for i, ch in enumerate(name):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        elif depth == 0 and i > 0 and CUT.match(name, i):
            return name[:i].strip(" -,")
    return name.strip()


def protein_name(known: list[str]) -> str | None:
    """The name mpstruc most often gives this protein, from the raw names of
    its existing entries (most recent first wins a tie)."""
    stems = [stem(n) for n in known if n and n.strip()]
    if not stems:
        return None
    count = Counter(s.lower() for s in stems)
    best = max(count.values())
    return next(s for s in stems if count[s.lower()] == best)


def build(det: dict | None, core: set[str], related: list[tuple[str, dict | None]] = (),
          known_names: list[str] = ()) -> dict | None:
    """Draft mpstruc fields for one structure.

    `related` lists (pdb, details) of structures from the same study that are
    suggested as its related entries; they are added to the description as
    mpstruc does ('with bound X, 3.10 Å: 9ABC'). `known_names` are the raw
    names of mpstruc entries for the same protein, most recent first; the
    draft names the protein the way they do. Text fields are mpstruc raw
    markup (entities escaped, <br />, <sub> and <sup> allowed).
    """
    if not det:
        return None
    ents = [e for e in det.get("polymer_entities") or []
            if (e.get("entity_poly") or {}).get("rcsb_entity_polymer_type") == "Protein"]
    mine = [e for e in ents if set((e.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []) & core]
    mine = mine or [e for e in ents if not ANTIBODY.search((e.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "")] or ents
    names = []
    for e in mine:
        d = _clean_name(((e.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip())
        if d and d not in names:
            names.append(d)
    orgs = [o for e in mine for o in e.get("rcsb_entity_source_organism") or []]
    hosts = [o.get("ncbi_scientific_name") for e in mine for o in e.get("rcsb_entity_host_organism") or []]
    species = next((o.get("ncbi_scientific_name") for o in orgs if o.get("ncbi_scientific_name")), None)
    lineage = {l.get("name") for o in orgs for l in o.get("taxonomy_lineage") or []}
    domain = next((d for d in DOMAINS if d in lineage), None)
    res = ((det.get("rcsb_entry_info") or {}).get("resolution_combined") or [None])[0]
    esc = lambda s: re.sub(r"&lt;(/?)(sub|sup)&gt;", r"<\1\2>", html.escape(s, quote=False))
    lines = []
    method = _method([m.get("method", "") for m in det.get("exptl") or []])
    if method:
        lines.append(method)
    pt = parts(det, core)
    for pdb, rd in related:
        rp = parts(rd, core) if rd else None
        line = esc(phrase(rp, pt)) if rp else ""
        rres = ((rd or {}).get("rcsb_entry_info") or {}).get("resolution_combined") or [None]
        res_txt = f"{rres[0]:.2f} Å" if rres[0] else ""
        lines.append(", ".join(x for x in (line, res_txt) if x) + f": {pdb}")
    known = protein_name(list(known_names))
    base = known or (esc(_cap(" / ".join(names[:2]))) if names else None)
    cond = esc(phrase(pt))
    name = (f"{base} {cond}" if cond else base) if base else None
    cit = det.get("rcsb_primary_citation") or {}
    pm = cit.get("pdbx_database_id_PubMed")
    pages = "-".join(dict.fromkeys(p for p in (cit.get("page_first"), cit.get("page_last")) if p)) or None
    return {
        "name": name,
        "species": esc(species) if species else None,
        "taxonomicDomain": domain,
        "expressedInSpecies": esc(hosts[0]) if hosts and hosts[0] else None,
        "resolution": f"{res:.2f}".rstrip("0").rstrip(".") if res else None,
        "description": "<br /> ".join(lines) or None,
        "bibliography": {
            "pubMedId": str(pm) if pm and int(pm) > 0 else None,
            "authors": esc(_authors(cit.get("rcsb_authors")) or "") or None,
            "year": str(cit["year"]) if cit.get("year") else None,
            "title": esc(cit["title"]) if cit.get("title") else None,
            "journal": esc(cit["journal_abbrev"]) if cit.get("journal_abbrev") and cit["journal_abbrev"] != "To Be Published" else None,
            "volume": cit.get("journal_volume"), "issue": cit.get("journal_issue"), "pages": pages,
            "doi": cit.get("pdbx_database_id_DOI"),
        },
        "pdbTitle": (det.get("struct") or {}).get("title"),
    }
