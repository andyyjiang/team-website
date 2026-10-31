"""Draft mpstruc entries for new PDB structures.

Builds the fields of an mpstruc <protein>/<memberProtein> record (name,
species, taxonomic domain, expression host, resolution, description and
bibliography) from the PDB, written the way mpstruc writes them:
"cryo-EM structure" on the first description line, "in complex with ..." and
"with bound ..." phrases, related structures as "..., 3.1 Å: 9ABC", and
authors as "Ahmed S, Huang M, & Lee CH". Curators should still edit the
name, which the PDB words differently from mpstruc.
"""
from __future__ import annotations

import html
import re

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


def _ligands(det: dict) -> list[str]:
    out = []
    for n in det.get("nonpolymer_entities") or []:
        comp = ((n.get("nonpolymer_comp") or {}).get("chem_comp") or {})
        cid, name = comp.get("id") or "", comp.get("name") or ""
        if not cid or cid in ADDITIVES:
            continue
        if cid in IONS:
            label = IONS[cid]
        elif len(name) <= 24 and not re.search(r"\d,\d|\(\d", name):
            label = name.lower() if name.isupper() else name
        else:
            label = cid
        if label not in out:
            out.append(label)
    return out[:4]


def _join(items: list[str]) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def condition(det: dict, core: set[str]) -> str:
    """The short phrase mpstruc uses for a structure: 'in complex with X, with bound Y'."""
    partners = []
    for e in det.get("polymer_entities") or []:
        if (e.get("entity_poly") or {}).get("rcsb_entity_polymer_type") != "Protein":
            continue
        ups = set((e.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or [])
        desc = ((e.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()
        if desc and not (ups & core) and desc not in partners:
            partners.append(desc)
    gprot = [p for p in partners if re.match(r"guanine nucleotide-binding protein", p, re.I)]
    if len(gprot) >= 2:
        partners = [p for p in partners if p not in gprot]
        partners.insert(0, "G protein heterotrimer")
    elif gprot:
        partners = [("G protein " + re.sub(r".*subunit\s*", "", gprot[0]).strip() if p in gprot else p) for p in partners]
    bits = []
    if partners:
        bits.append("in complex with " + _join(partners[:3]))
    ligs = _ligands(det)
    if ligs:
        bits.append("with bound " + _join(ligs))
    return ", ".join(bits)


def build(det: dict | None, core: set[str], related: list[tuple[str, dict | None]] = ()) -> dict | None:
    """Draft mpstruc fields for one structure.

    `related` lists (pdb, details) of structures from the same study that are
    suggested as its related entries; they are added to the description as
    mpstruc does ('with bound X, 3.1 Å: 9ABC'). Text fields are mpstruc raw
    markup (entities escaped, <br /> and <sup> allowed).
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
    esc = lambda s: html.escape(s, quote=False)
    lines = []
    method = _method([m.get("method", "") for m in det.get("exptl") or []])
    if method:
        lines.append(method)
    cond = condition(det, core)
    if cond:
        lines.append(esc(cond).replace("&lt;sup&gt;", "<sup>").replace("&lt;/sup&gt;", "</sup>"))
    for pdb, rd in related:
        rcond = condition(rd, core) if rd else ""
        rres = ((rd or {}).get("rcsb_entry_info") or {}).get("resolution_combined") or [None]
        phrase = esc(rcond).replace("&lt;sup&gt;", "<sup>").replace("&lt;/sup&gt;", "</sup>") or "structure"
        lines.append(f"{phrase}{f', {rres[0]:.2f} &Aring;' if rres[0] else ''}: {pdb}")
    cit = det.get("rcsb_primary_citation") or {}
    pm = cit.get("pdbx_database_id_PubMed")
    pages = "-".join(dict.fromkeys(p for p in (cit.get("page_first"), cit.get("page_last")) if p)) or None
    return {
        "name": esc(_cap(" / ".join(names[:2]))) if names else None,
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
