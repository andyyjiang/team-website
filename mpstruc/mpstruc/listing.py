"""Draft mpstruc entries for new PDB structures.

Builds the fields of an mpstruc <protein>/<memberProtein> record (name,
species, taxonomic domain, expression host, resolution, description and
bibliography) from the PDB, written the way mpstruc curators write them.
The conventions below were measured on mpstruc's own text (2024-2026
entries of the 2 and 15 Oct 2026 snapshots):

* The NAME carries the structure's condition: "<protein>, apo form",
  "<protein> - G<sub>s</sub> complex with bound X", "<protein> in complex
  with Y". For a daughter, <protein> is the name mpstruc already uses for
  that protein (its parent's family), not the PDB's entity name.
* The DESCRIPTION's first line is only the method ("cryo-EM structure");
  each related structure follows on its own line as
  "<what differs>, 3.10 &Aring;: 9ABC" (two decimals in 99% of lines).
* A related line names only what DIFFERS within the study: a ligand, partner,
  mutation, pH, state or map that every structure of the group shares is
  left out. Facets are written in the order curators use: mutation,
  complex, bound ligand, environment, pH, temperature, state, class/map.
  Two structures of a study never get the same line when their PDB titles
  differ: the title words that set them apart are added.
* Names and descriptions get mpstruc's <sub>/<sup> markup (markup.py).
  Bibliography titles stay as the journal printed them: none of the 803
  titles curated for 2024-2026 papers carries markup.

Credit: parts of this come from PR #4 (branch mpstruc-listing-language, a
collaborator's session), each kept only because it improved the measured
match with curators' text: cutting a protein name at 'in a ...' / 'in DMPC
...' and keeping entities whole (its stem()); a ligand with only a systematic
name taking the name the PDB title gives it; reusing the name of mpstruc
entries with exactly the same proteins when the family gives none (its
known_names, in automation.add_listings). The three gaps its description
listed are fixed here too: the most common spelling of a family name wins,
BRIL / soluble cytochrome b562 and anti-BRIL or anti-Fab antibodies are never
partners, and detergents such as CHAPSO are never ligands. Its 40-residue
peptide cut-off was tried and not kept (100 matches curators better).
"""
from __future__ import annotations

import copy
import html
import json
import re
import time
from collections import Counter

from .markup import format_text

DOMAINS = ("Eukaryota", "Bacteria", "Archaea", "Viruses")
# Ions written the way mpstruc writes them, and the words a PDB title uses for them.
IONS = {"CA": "Ca<sup>2+</sup>", "MG": "Mg<sup>2+</sup>", "ZN": "Zn<sup>2+</sup>", "MN": "Mn<sup>2+</sup>",
        "NA": "Na<sup>+</sup>", "K": "K<sup>+</sup>", "CL": "Cl<sup>-</sup>", "CU": "Cu<sup>2+</sup>",
        "CU1": "Cu<sup>+</sup>", "FE": "Fe<sup>3+</sup>", "FE2": "Fe<sup>2+</sup>", "CD": "Cd<sup>2+</sup>",
        "LI": "Li<sup>+</sup>", "RB": "Rb<sup>+</sup>", "CS": "Cs<sup>+</sup>", "BA": "Ba<sup>2+</sup>",
        "SR": "Sr<sup>2+</sup>", "NI": "Ni<sup>2+</sup>", "CO": "Co<sup>2+</sup>", "TL": "Tl<sup>+</sup>",
        "GD": "Gd<sup>3+</sup>", "YB": "Yb<sup>3+</sup>", "HG": "Hg<sup>2+</sup>", "PT": "Pt<sup>2+</sup>"}
ION_WORDS = {"CA": r"calcium|ca2\+|ca\b", "MG": r"magnesium|mg2\+", "ZN": r"zinc|zn2\+", "MN": r"manganese|mn2\+",
             "NA": r"sodium|na\+", "K": r"potassium|k\+", "CL": r"chloride|cl-", "CU": r"copper|cu2\+",
             "CU1": r"copper|cu\+", "FE": r"iron|fe3\+", "FE2": r"iron|fe2\+", "CD": r"cadmium|cd2\+",
             "LI": r"lithium|li\+", "RB": r"rubidium|rb\+", "CS": r"caesium|cesium|cs\+", "BA": r"barium|ba2\+",
             "SR": r"strontium|sr2\+", "NI": r"nickel|ni2\+", "CO": r"cobalt|co2\+", "TL": r"thallium|tl\+",
             "GD": r"gadolinium|gd3\+", "YB": r"ytterbium|yb3\+", "HG": r"mercury|hg2\+", "PT": r"platinum|pt2\+"}
# 'Potassium channel', 'sodium-dependent transporter': an ion word naming the protein, not a bound ion.
ION_NOUN = (r"(?:channels?|transporters?|pumps?|ATPases?|exchangers?|symporters?|antiporters?|cotransporters?|"
            r"uniporters?|activated|dependent|gated|selective|sensing|transporting|release|permeable|coupled|"
            r"sensitive|leak)")
# Ligands mpstruc writes by a short common name.
COMMON = {"ATP": "ATP", "ADP": "ADP", "AMP": "AMP", "ANP": "AMP-PNP", "ACP": "AMPPCP", "AGS": "ATP&gamma;S",
          "GTP": "GTP", "GDP": "GDP", "GSP": "GTP&gamma;S", "GNP": "GMPPNP", "BEF": "BeF<sub>3</sub><sup>-</sup>",
          "ALF": "AlF<sub>4</sub><sup>-</sup>", "MGF": "MgF<sub>3</sub><sup>-</sup>", "VO4": "vanadate",
          "SRO": "serotonin", "GLU": "glutamate", "GLY": "glycine", "HSM": "histamine", "ACH": "acetylcholine",
          "EPE": "HEPES", "CMP": "cAMP", "PCG": "cGMP", "I3P": "IP<sub>3</sub>", "PIO": "PIP<sub>2</sub>",
          "CLR": "cholesterol", "Y01": "CHS", "RET": "retinal", "HEM": "heme", "SAM": "SAM", "SAH": "SAH",
          "NAD": "NAD<sup>+</sup>", "NAI": "NADH", "FAD": "FAD", "FMN": "FMN", "UQ": "ubiquinone",
          "8PE": "phosphatidylethanolamine", "TTX": "tetrodotoxin", "COA": "CoA", "ACO": "acetyl-CoA",
          "SUC": "succinate", "CIT": "citrate", "GLC": "glucose", "ASP": "aspartate"}
# Amino acids as bound ligands: mpstruc writes the three-letter code ('with bound Leu').
COMMON.update({a.upper(): a for a in "Ala Arg Asn Cys Gln His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val".split()})
# Water, buffers, cryoprotectants, detergents and common lipids: never "bound ligands" in mpstruc
# unless the PDB title names them.
ADDITIVES = set("""HOH DOD GOL EDO PEG PG4 PGE 1PE P6G PE4 PE5 SO4 PO4 ACT FMT DMS MPD TRS EPE MES BME IMD
CIT TLA ACY IPA EOH MOH NO3 IOD BR UNX UNL BOG LMT LMN LDA C8E BNG HTG OLC OLA OLB 78M 78N 7E9 MPG
Y01 CLR PCW POV PEE PLC LPP PGW CDL 3PE PC1 6PL P5S PEF PEV D10 D12 OCT HEX DD9 R16 TRD LFA PLM MYR
STE PA8 PA6 PTY NAG NDG MAN BMA FUC GAL GLC SIA DMU UMQ CPS SDS DPC LDM 9Z9 AV0 PX4 PCF PGT U10
7PH PTD 8PE 44E 2CV 17F J40 P0E LHG LMG DGD SQD AJP LBN 6OU PEK PSF CHD A85 PIE EPH 3TL NA1 POP
SPD SPM PUT DTT TCE EDT EGT AZI SCN FLC MLI MLA PGO PG0 BU3 2PE 12P 15P PE3 XPE MGE LNL DGG LMU
1N7 C15 TRT CXE CE1 CE9 SOG BGL LMO RWB P4K
""".split())
# Cofactors and pigments that every structure of a complex carries: named only when the title names them.
COFACTORS = set("""HEM HEC HEA HEB CLA CHL BCL BPH BCR CRT LUT NEX XAT PQN PL9 FAD FMN NAD NAP NDP NAI SF4
FES F3S MQ7 MQ8 UQ8 UQ1 UQ2 UQ6 UQ9 DGA OEX BCB A86 KC1 ECH LMU CU1 RET FE2 MQ9 3PH WRA CHL ZEX
PHO BCT CUA CU HEO HAS 3PE SMA DMV MQ9 CYC PEB PUB BLA FDA FAD CMO
PER PEO OXY""".split())  # peroxide, H2O2, O2: oxidase active-site ligands
# Other spellings of common ligands, so the title's own spelling is used ('AMPPNP' vs 'AMP-PNP').
COMMON_ALT = {"ANP": ["AMPPNP", "AMP-PNP", "AppNHp", "ANP"], "ACP": ["AMPPCP", "AMP-PCP", "AppCp"],
              "AGS": ["ATPgammaS", "ATP-gamma-S", "ATPγS"], "GSP": ["GTPgammaS", "GTPγS"],
              "GNP": ["GppNHp", "GMPPNP"], "BEF": ["BeFx", "BeF3", "BeF3-"], "ALF": ["AlFx", "AlF4", "AlF4-"],
              "SRO": ["5-HT", "serotonin"], "I3P": ["IP3", "InsP3"], "PIO": ["PIP2", "PI(4,5)P2"],
              "CMP": ["cAMP"], "PCG": ["cGMP"], "GLU": ["glutamate", "L-glutamate", "Glu"]}
LIPIDISH = re.compile(r"\b\d{1,2}:\d\b|\b(?:lpa|lpc|lpe|lps|pip\d?|dag)\b(?=.*\d)|(?:hex|hept|oct|non|dec|undec|dodec|tridec|tetradec|pentadec|hexadec|heptadec|octadec)ane\b|"
                      r"phosphatidyl|phospho(?:choline|ethanolamine|serine|glycerol|inositol)|glycero-?3|sn-glycer|oleoyl|palmitoyl|stearoyl|myristoyl|cardiolipin|maltoside|glucoside|"
                      r"digitonin|hemisuccinate|lauryl|dodecyl|undecyl|decyl|nonyl|octyl|neopentyl|cholest|"
                      r"diacylglycerol|lipid|detergent|sterol|polyethylene|glycol|galactosyl|"
                      # zwitterionic and other detergents (PR #4 listed CHAPSO drafted as a ligand)
                      r"\bchapso?\b|cholamidopropyl|ammonio-1-propanesulfonate|fos-?choline|\btriton\b|\btween\b|zwittergent|"
                      r"cymal|\bdeoxycholate\b|\bcholate\b|"
                      # fatty acids: 8 carbons or more ('hexacosanoic', 'alpha-linolenic'), not 'ethanoic'
                      r"\b[a-z-]*(?:octan|nonan|decan|cosan|octadec|hexadec|tetradec|dodec|icos|docos)[a-z]*oic acid\b|"
                      r"\b[a-z-]*(?:oleic|stearic|palmitic|myristic|lauric|linoleic|linolenic|arachidonic|petroselinic|"
                      r"elaidic|vaccenic) acid\b", re.I)
ANTIBODY = re.compile(r"\b(fab|fabs|nanobody|nanobodies|nb\d*|scfv\d*|antibody|heavy chain|light chain|sybody|"
                      r"darpin|megabody|igg\d*|vhh|monobody|fv|fragment antigen|antigen[- ]binding|"
                      r"single[- ]chain variable|legobody|anti-?fab)\b", re.I)
STOP = set("a an and the of in on at to by for from with without into via as is under its their".split())


# ---------------------------------------------------------------- small helpers

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


def _esc(s: str) -> str:
    """Escape text but keep the <sub>/<sup> tags and entities this module writes itself.

    Only a properly closed pair ('Ca<sup>2+</sup>') becomes markup again; a stray
    or unclosed tag from PDB text is dropped, so it can never subscript the rest
    of a card or leave the XML download unbalanced."""
    s = html.escape(s, quote=False)
    s = re.sub(r"&amp;([A-Za-z]+);", r"&\1;", s)
    s = re.sub(r"&lt;(sub|sup)&gt;((?:[^&<]|&[A-Za-z]+;){1,20}?)&lt;/\1&gt;", r"<\1>\2</\1>", s)
    return re.sub(r"&lt;/?(?:sub|sup)&gt;", "", s)


def _mp(s: str) -> str:
    """Escaped text with mpstruc's sub/superscripts added."""
    return format_text(_esc(s)) or ""


def _join(items: list[str]) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _norm(s: str) -> str:
    """Lowercase, markup-free form for comparisons."""
    s = html.unescape(re.sub(r"<[^>]+>", "", s or ""))
    return re.sub(r"[^a-z0-9α-ω+]+", " ", s.lower()).strip()


FUSION_NAME = re.compile(r"^(green fluorescent protein|gfp|egfp|t4 lysozyme|endolysin|lysozyme|soluble cytochrome b562|"
                         r"bril|maltose(?:/maltodextrin)?[- ]binding.*|maltodextrin[- ]binding.*|thioredoxin.*|"
                         r"rubredoxin|glycogen synthase.*|mcherry|venus|fusion protein|(?:engineered )?.*luciferin.*|"
                         r".*luciferase.*|(?:large|small)? ?bit|nanobit|lgbit|smbit|hibit|nanoluc)$", re.I)


# Fusion tags and crystallisation chaperones (as automation.FUSION_UNIPROT).
FUSION_TAG_UNIPROT = {"P00720", "P0ABE7", "P0AEX9", "P42212", "P0AA25", "Q9ZEU2", "P61626", "P63165", "P0ACJ8",
                      "P02185"}
# Of those, the ones that are never a real binding partner when they come as a chain of their
# own. Maltose-binding protein (P0AEX9) is the substrate-binding partner of the maltose
# transporter ('MalFGK<sub>2</sub>-MBP'), and thioredoxin (P0AA25) has partners of its own.
TAG_ONLY_UNIPROT = FUSION_TAG_UNIPROT - {"P0AEX9", "P0AA25"}
FUSION_TAG_GENES = {"gfp", "cybc", "male", "e", "trxa", "rubr", "lyz", "sumo1", "mb"}
# BRIL (soluble cytochrome b562) written into a fused protein's name: 'OaVMAT2-BRIL',
# 'Adenosine receptor A2a/Soluble cytochrome b562 chimera', 'NorA-BRIL(3A) fusion'.
BRIL_IN_NAME = re.compile(r"\s*(?:[-/,]|\band\b)?\s*(?:E\.? ?coli )?(?:soluble |apo)?cytochrome b562(?:RIL)?"
                          r"(?:\s+(?:complex|chimera|fusion(?: protein)?))?|(?<!anti)(?<!anti-)(?<!anti )-?\bBRIL\b(?:\(\w+\))?(?:\s+fusion)?", re.I)
# Antibodies that only hold the fusion tag or another antibody (anti-BRIL Fab BAG2, anti-Fab nanobody).
HELPER_ANTIBODY = re.compile(r"anti-?\s?bril|anti-?\s?fab|\bBAG2\b|SRP2070|scfv[-\s]?16|\bnb[-\s]?35\b|fab[-\s]?(?:cb6|16)\b",
                             re.I)


def _is_fusion_tag(e: dict) -> bool:
    """A chain that is only a fusion tag or chaperone (a separate BRIL chain),
    never a partner. Not maltose-binding protein or thioredoxin, which can be
    real partners, and not a chain that fuses a tag to something else
    ('Maltodextrin-binding protein,Autotransporter EspP chimera')."""
    desc, ups = _desc(e), _uniprots(e)
    if re.search(r"malto|thioredoxin", desc, re.I) or (ups and not ups <= TAG_ONLY_UNIPROT):
        return False
    parts = [re.sub(r"\s+(?:chimera|fusion(?: protein)?)$", "", x, flags=re.I).strip() for x in desc.split(",")]
    parts = [x for x in parts if x]
    return bool(parts) and all(FUSION_NAME.match(x) or BRIL_IN_NAME.fullmatch(x) for x in parts)


def _clean_name(desc: str) -> str:
    """PDB entity name -> mpstruc-like protein name: drop fusion partners
    ('X,Green fluorescent protein', 'X/Soluble cytochrome b562 chimera',
    'X-BRIL') and 'Isoform 1 of'."""
    if BRIL_IN_NAME.search(desc) and not FUSION_NAME.match(desc.strip()):
        bare = re.sub(r"^and\s+|\s+and$", "", BRIL_IN_NAME.sub("", desc).strip(" ,/-"), flags=re.I).strip(" ,/-")
        # Each '/'- or ','-separated part loses its own 'chimera' or 'fusion', and a part
        # repeated around the tag goes ('A2a/b562/A2a chimera,b562,A2a' -> 'A2a').
        safe = re.sub(r"(?<![A-Za-z0-9])\d+(?:,\d+)+(?=[-\s)'])", lambda m: m[0].replace(",", "\0"), bare)
        pieces, seen = re.split(r"(\s*[/,]\s*)", safe), set()
        bare = ""
        for i in range(0, len(pieces), 2):
            x = re.sub(r"^(?:chimera|fusion) protein of (?:human )?|\s+(?:chimera|fusion(?: protein)?)$", "",
                       pieces[i].strip(), flags=re.I).strip()
            if x and x.lower() not in seen:
                seen.add(x.lower())
                bare += (pieces[i - 1].strip() if bare else "") + x
        bare = bare.replace("\0", ",")
        if len(bare) >= 3:
            desc = bare
    rep = re.fullmatch(r"(.+?)(?:,\s*\1)+(,.*)?", desc)  # 'X 1,X 1,Fusion protein' -> 'X 1,Fusion protein'
    if rep:
        desc = rep[1] + (rep[2] or "")
    # Commas of a locant ('Inositol 1,4,5-trisphosphate', 'beta-(1,2)-glucan') do not separate parts.
    safe = re.sub(r"(?<![A-Za-z0-9])\d+(?:,\d+)+(?=[-\s)'])", lambda m: m[0].replace(",", "\0"), desc)
    parts = list(dict.fromkeys(p.strip().replace("\0", ",") for p in safe.split(",") if p.strip()))  # 'X,X,X' -> 'X'
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


def _res(det: dict | None) -> float | None:
    return (((det or {}).get("rcsb_entry_info") or {}).get("resolution_combined") or [None])[0]


def _title(det: dict | None) -> str:
    return ((det or {}).get("struct") or {}).get("title") or ""


def _protein_entities(det: dict) -> list[dict]:
    return [e for e in det.get("polymer_entities") or []
            if (e.get("entity_poly") or {}).get("rcsb_entity_polymer_type") == "Protein"]


def _uniprots(e: dict) -> set[str]:
    return set((e.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or [])


def _desc(e: dict) -> str:
    return ((e.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()


def _length(e: dict) -> int | None:
    return (e.get("entity_poly") or {}).get("rcsb_sample_sequence_length")


# ---------------------------------------------------------------- names the title gives

ROLE = re.compile(r"^(?:the |a |an )?(?:(?:partial|full|inverse|allosteric|positive|negative|biased|small[- ]molecule|"
                  r"peptide|endogenous|natural|selective|covalent|irreversible|novel)\s+)*"
                  r"(?:agonists?|antagonists?|inhibitors?|activators?|modulators?|blockers?|openers?|PAMs?|NAMs?|"
                  r"ago-PAMs?|substrates?|drugs?|ligands?|compounds?|molecules?)\s+", re.I)
TITLE_GENERIC = re.compile(
    r"^(?:proteins?|receptors?|(?:super|sub)?complex(?:es)?|nuclear|pores?|cores?|rings?|baskets?|bod(?:y|ies)|"
    r"basal|inner|outer|motors?|assembl(?:y|ies)|structures?|apo|lipids?|nanodiscs?|detergents?|ligands?|substrates?|"
    r"agonists?|antagonists?|inhibitors?|cryo-?em|cryoem|cryo-electron microscopy|x-ray|crystal|human|mouse|rat|"
    r"bovine|yeast|bacterial|full[- ]length|wild[- ]type|wt|mutant|state|form|conformation|the|its|it|an?|"
    r"g[- ]?proteins?|heterotrimer(?:ic)?(?: g[- ]?proteins?)?|g?gustducin|transducin|(?:mini-?)?g(?:s|i|o|q|t|olf|11|12|13)\d?|"
    r"(?:beta|β)-?arrestin-?\d?|arrestin|mini|nanobody\s*\d*|nb\d+|scfv\s*\d*|fab\s*\w*|antibody|calcium|sodium|"
    r"potassium|zinc|magnesium|chloride|ions?|water|membranes?|liposomes?|amphipols?|resolution|angstrom|"
    r"\d+(?:\.\d+)?\s*(?:a|å|angstrom)?|ph\s*[\d.]+)$", re.I)


PROTEINISH = re.compile(r"\b(?:receptors?|channels?|transporters?|proteins?|enzymes?|pumps?|complex|subunits?|"
                        r"domains?|[a-z]+ases?|TMD|LBD|NTD|ECD|CTD|ICD|state|structure)$|^(?:TMD|LBD|NTD|ECD|CTD)\b", re.I)


def title_names(title: str, aliases: set[str] = frozenset()) -> list[str]:
    """Ligand and partner names a PDB title gives, in title order:
    'mGlu2 bound to agonist LY354740 and PAM JNJ-46281222' -> ['LY354740', 'JNJ-46281222'],
    'Human Slo1-Charybdotoxin complex' -> ['Charybdotoxin'] (Slo1 being the protein itself)."""
    t = re.sub(r"\s+", " ", title or "")
    found: list[str] = []
    stop = (r"(?=\s+(?:in|at|under|obtained|from|determined|using|reveal\w*|via|bound|state|structure|with|"
            r"complex|complexes|"
            r"conformation|form|mode|site|pocket|reconstituted|solubili[sz]ed|embedded|captured|after|before|"
            r"during|without)\b|[,;()]|\s+-\s+|$)")
    for m in re.finditer(r"\b(?:bound (?:to|with)|in complex with|complexed with|complex with|with|activated by|"
                         r"inhibited by|occupied by|and|complex structure of|in (?:the )?presence of)\s+(.+?)" + stop,
                         t, re.I):
        for item in re.split(r"\s+and\s+|\s*/\s*|,\s*|\s+&\s+", m[1]):
            found.append(item)
    for m in re.finditer(r"([\w.+]+(?:-[\w.+]+)+)(?=\s+(?:protein\s+)?complex\b)", t):
        found.append(m[1])
    for m in re.finditer(r"\b([\w.+-]+)-(?:bound|occupied|loaded)\b", t):
        found.append(m[1])
    m = re.match(r"([A-Za-z0-9][\w.+-]*)\s+(?:bound|occupied)\b", t)  # 'AR234958 bound Mas1 receptor'
    if m:
        found.append(m[1])
    # 'Slo1-paxilline', 'MCT2-embigin-AR-C155858': a hyphen chain that includes the
    # protein itself or its G protein lists the other members.
    split: list[str] = []
    for item in found:
        parts = _chain(item)
        if len(parts) > 1 and any(_is_alias(p, aliases) or TITLE_GENERIC.match(p) for p in parts):
            split += parts
        else:
            split.append(item)
    out: list[str] = []
    for item in split:
        item = ROLE.sub("", re.sub(r"^(?:(?:the|a|an|human|mouse|rat|bovine|bound)\s+)+", "", item.strip(" ."),
                                   flags=re.I)).strip(" .")
        if not item or len(item) < 2 or len(item.split()) > 4 or len(item) > 40 or TITLE_GENERIC.match(item):
            continue
        if MUTATION.fullmatch(item) or re.fullmatch(r"\d+|[a-z]+", item) and item.lower() in STOP:
            continue
        if _is_alias(item, aliases) or item in out or PROTEINISH.search(item):
            continue
        out.append(item)
    return out


def _chain(item: str) -> list[str]:
    """'MCT2-embigin-AR-C155858' -> ['MCT2', 'embigin', 'AR-C155858'] (a 1-2 capital
    prefix such as 'AR-' or 'L-' belongs to the next piece, as does a bare number)."""
    parts: list[str] = []
    for piece in item.split("-"):
        if parts and (re.fullmatch(r"[A-Z]{1,2}|\d+", parts[-1]) or re.fullmatch(r"\d+", piece)):
            parts[-1] += "-" + piece
        else:
            parts.append(piece)
    return [p for p in parts if p]


def _is_alias(item: str, aliases: set[str]) -> bool:
    """Is `item` just another name of the structure's own protein ('Slo1' for
    the gene SLO, 'NET' for '... transporter (NET)')? A word that merely starts
    the protein's name ('C5a' in 'C5a receptor') is not."""
    k = _norm(item).replace(" ", "")
    if not k:
        return True
    for a in aliases:
        if re.match(rf"{re.escape(_norm(item))} (?:type \d+ )?receptor", _norm(a)):
            continue  # 'neurotensin' is the ligand of 'Neurotensin receptor', not the receptor
        a = _norm(a).replace(" ", "")
        if len(a) < 2:
            continue
        if k == a or k.rstrip("r") == a or (len(a) >= 3 and k.startswith(a) and len(k) - len(a) <= 2) \
                or (len(k) >= 4 and a.startswith(k.rstrip("r"))):
            return True
    return False


# ---------------------------------------------------------------- partners

def _g_label(descs: list[str]) -> str:
    """The G protein of a complex as mpstruc writes it ('G<sub>s</sub>', 'miniG<sub>q</sub>')."""
    alpha = [d for d in descs if not re.search(r"beta|gamma", d, re.I)]
    if not alpha:
        return "G protein"
    for d in alpha:
        m = re.search(r"mini[-\s]?G\(?(s|q|i|o|12|13)\)?(?![a-z])", d, re.I)
        if m:
            return f"miniG<sub>{m[1].lower()}</sub>"
    for d in alpha:
        m = (re.search(r"G\((s|i|o|q|t|olf|z|11|12|13)\)", d, re.I)
             or re.search(r"\bG[-\s]?alpha[-\s]?(s|i|o|q|t|olf|z|11|12|13)", d, re.I)
             or re.search(r"\bGNA(S|I|O|Q|T)", d)
             or re.search(r"subunit alpha-(1[1-6])\b", d)
             or re.search(r"(?<![A-Za-z0-9])G(s|i|o|q|t|olf)\d?(?![A-Za-z0-9])", d, re.I))
        if m:
            return f"G<sub>{m[1].lower()}</sub>"
    return "G protein"


def _title_g(title: str) -> str | None:
    """A G protein or arrestin partner the title names but the model may lack."""
    m = re.search(r"(?:beta|β)[- ]?arrestin[- ]?(\d)", title, re.I)
    if m:
        return f"β-arrestin {m[1]}"
    if re.search(r"\barrestin\b", title, re.I):
        return "arrestin"
    if re.search(r"gustducin", title, re.I):  # mpstruc: 'G<sub>g</sub> complex'
        return "G<sub>g</sub>"
    if re.search(r"transducin", title, re.I):
        return "G<sub>t</sub>"
    m = re.search(r"mini[-\s]?G(s/q|s/i|s|q|i|o|12|13)(?![a-z])", title, re.I)
    if m:
        return f"miniG<sub>{m[1].lower()}</sub>"
    m = re.search(r"(?<![A-Za-z0-9])G(s|i|o|q|t|olf|z)(\d?)(?![A-Za-z0-9])", title)
    if m and re.search(r"complex|coupled|bound|protein|-G", title, re.I):
        return f"G<sub>{m[1].lower()}{m[2]}</sub>"
    if re.search(r"\bG[- ]?proteins?\b(?![- ]coupled)|heterotrimeric G", title, re.I):
        return "G protein"
    return None


def _partner_label(desc: str) -> str:
    d = _clean_name(desc)
    d = re.sub(r",? (?:mitochondrial|chloroplastic|isoform \w+)$", "", d, flags=re.I)
    if re.fullmatch(r"Calmodulin(?:-\d)?", d, re.I):
        return "calmodulin"
    # 'Embigin' -> 'embigin', but keep acronyms and gene-style names ('TARP', 'Kv1.2').
    if re.fullmatch(r"[A-Z][a-z]{3,}(?: [a-z][a-z-]*)*", d):
        return d[0].lower() + d[1:]
    return d


def _is_antibody(desc: str) -> bool:
    return bool(ANTIBODY.search(desc))


def _genes(e: dict) -> list[str]:
    return [g["value"] for o in e.get("rcsb_entity_source_organism") or [] for g in o.get("rcsb_gene_name") or []
            if g.get("value")]


GENERIC_WORD = re.compile(r"^(?:protein|proteins|subunit|chain|domain|family|member|receptor|channel|toxin|complement|"
                          r"alpha|beta|gamma|delta|component|containing|binding|type|like|putative|"
                          r"uncharacteri[sz]ed|precursor|fragment|peptide|isoform|human|potassium|sodium|calcium|"
                          r"complex|complexes|supercomplex|subcomplex|nuclear|pore|core|ring|basket|body|basal|inner|"
                          r"outer|engineered|terminal|n-terminal|c-terminal|structure|structural|regulatory|"
                          r"catalytic|accessory|associated|factor|unit|large|small|light|heavy|mitochondrial|"
                          r"chloroplastic|transporter|system|assembly|motor|state|form|bound|with)$", re.I)


def _said(e: dict, title: str, own: str = "") -> str | None:
    """How the title names a partner chain (its gene symbol or a distinctive
    word of its name that is not also a word of the protein's own name)."""
    own_words = set(_norm(own.replace("|", " ")).split())
    words = [w for w in re.findall(r"[A-Za-z0-9][\w-]*", _desc(e))
             if len(w) >= 3 and w.lower() not in STOP and not GENERIC_WORD.match(w) and _norm(w) not in own_words]
    # Distinctive tokens first ('CX3CL1', 'Nup155'), then plain words ('fractalkine').
    words.sort(key=lambda w: not re.search(r"\d|.[A-Z]", w))
    return _in_title(_genes(e) + words, title)


def partners(det: dict, core: set[str], title: str | None = None, own: str = "") -> dict:
    """Non-core protein chains: the G protein or arrestin, other partners,
    antibodies and peptide ligands (chains of 100 residues or fewer, which
    mpstruc writes as 'with bound X'). Each partner is (label, named in title)."""
    title = _title(det) if title is None else title
    out = {"g": None, "g_from": None, "partners": [], "antibodies": [], "peptides": []}
    g_descs = []
    for e in _protein_entities(det):
        ups, desc = _uniprots(e), _desc(e)
        if not desc or (ups & core) or _is_fusion_tag(e):
            continue  # (a separate BRIL chain is a tag, not a partner: a known gap of PR #4)
        if re.search(r"guanine nucleotide-binding protein|\bG(?:[soqitz]|olf|1[1-6])?[-\s]?proteins?\b(?![-\s]coupled)|\bG[-\s]?alpha|"
                     r"mini[-\s]?G[soqi]\b|\bG\(?[soqi]\)? (?:alpha|subunit)|\bGNA[ISOQT]|gustducin|transducin",
                     desc, re.I):
            g_descs.append(desc)
            continue
        m = re.match(r"(?:beta|β)[- ]arrestin[- ]?(\d)", _clean_name(desc), re.I)
        if m:
            out["g"] = f"β-arrestin {m[1]}"
            continue
        if _is_antibody(desc):
            m = re.search(r"\b(?:Nb|Sb|Fab|scFv|Mb)[-\s]?\w*\d\w*", desc, re.I)
            said = m and _in_title([m[0], m[0].replace(" ", "")], title)
            if HELPER_ANTIBODY.search(desc):
                # The standard G-protein helpers (Nb35, scFv16) and the anti-BRIL Fab or
                # anti-Fab nanobody say nothing about the structure (PR #4).
                said = None
            out["antibodies"].append((said or desc, bool(said)))
            continue
        tarp = re.match(r"voltage-dependent calcium channel gamma-([2-8]) subunit", desc, re.I)
        if tarp:  # mpstruc: 'TARP-γ2'
            said, lab = _in_title(["TARP"], title), f"TARP-γ{tarp[1]}"
        elif SEQUENCE.fullmatch(desc.upper()):
            said, lab = None, "peptide"  # 'LYS-LYS-PRO-TYR-...': a bare sequence, never shown as such
        elif _same_complex(desc, own):
            continue  # 'ABC transporter ATP-binding protein' of an ABC transporter: part of the protein itself
        else:
            said = _said(e, title, own)
            lab = _lower_name(said) if said else _partner_label(desc)
            genes = _genes(e)
            if not said and len(lab.split()) >= 4 and genes and re.fullmatch(r"[A-Z][A-Z0-9-]{2,9}", genes[0]):
                lab = genes[0]  # 'Potassium voltage-gated channel subfamily KQT member 3' -> 'KCNQ3'
            elif not said and len(lab.split()) >= 5:
                # 'N-terminal fragment of the engineered fractalkine CX3CL1.44' -> 'CX3CL1.44'
                tok = next((w for w in re.findall(r"[A-Za-z0-9][\w.-]*[A-Za-z0-9]", desc)
                            if len(w) >= 3 and re.search(r"\d|.[A-Z]", w) and not GENERIC_WORD.match(w)
                            and not re.fullmatch(r"[\d.]+|[NC]-terminal|k?Da", w, re.I)), None)
                lab = tok or lab
            if _norm(lab) in {_norm(a) for a in own.split("|")} or _norm(desc) in {_norm(a) for a in own.split("|")}:
                continue  # another copy of the protein itself
        kind = "peptides" if (_length(e) or 999) <= 100 else "partners"
        if lab not in (x for x, _ in out[kind]):
            out[kind].append((lab, bool(said)))
    tg = _title_g(title)
    if g_descs and not out["g"]:
        # The title's own subtype ('Gi3', 'miniGs/q') when it gives one.
        use_title = bool(tg and tg != "G protein" and "arrestin" not in tg)
        out["g"], out["g_from"] = (tg, "title") if use_title else (_g_label(g_descs), "chains")
    if not out["g"]:
        out["g"], out["g_from"] = tg, "title"
    return out


# A peptide named only by its residues ('LYS-LYS-PRO-TYR-ILE-LEU', 'ASP-LE1-PHE-...').
SEQUENCE = re.compile(r"(?:[A-Z0-9]{3}[- ]){3,}[A-Z0-9]{3}")
SUBUNIT_WORDS = re.compile(r"\b(?:atp|atpase|nucleotide|transmembrane|permease|protein|proteins|subunit|subunits|component|domain|"
                           r"periplasmic|substrate-binding|binding|chain|small|large|alpha|beta|gamma|delta|"
                           r"[a-z]|\d+)\b")


def _same_complex(desc: str, own: str) -> bool:
    """Is a partner chain another subunit of the protein itself? True when, without
    subunit words, its name ('ABC transporter') is part of one of the protein's
    own names ('... export ABC transporter permease subunit PmtD')."""
    head = re.sub(r"\s+", " ", SUBUNIT_WORDS.sub(" ", _norm(desc))).strip()
    if len(head.split()) < 2:
        return False
    for a in own.split("|"):
        a = re.sub(r"\s+", " ", SUBUNIT_WORDS.sub(" ", _norm(a))).strip()
        if a and f" {head} " in f" {a} ":
            return True
    return False


# ---------------------------------------------------------------- ligands

def _in_title(words: list[str], title: str) -> str | None:
    """The first of `words` the title uses, spelled as the title spells it."""
    for w in words:
        w = (w or "").strip()
        if len(w) < 2:
            continue
        m = re.search(rf"(?<![A-Za-z0-9]){re.escape(w)}(?![A-Za-z0-9])", title, re.I)
        if m:
            return m.group(0)
    return None


def _readable(name: str) -> bool:
    return 0 < len(name) <= 30 and not re.search(r"\d,\d|\(\d|\[|\{|~|\d[A-Z]?-[a-z]{2,}.*\d|_\d|\d{6,}|\*", name)


def ligands(det: dict, title: str | None = None) -> list[tuple[str, str, str]]:
    """(chem id, mpstruc label, how) for each bound ligand worth naming;
    `how` is 'title' (the title names it), 'name' (a readable PDB name) or
    'id' (only the chemical ID is known).

    Drugs and substrates are kept; ions, water, buffers, detergents, lipids
    and cofactors only when the PDB title names them. A ligand the title
    names is written the way the title writes it ('LSD', 'paxilline')."""
    title = title if title is not None else _title(det)
    out: list[tuple[str, str, str]] = []
    for n in det.get("nonpolymer_entities") or []:
        comp = ((n.get("nonpolymer_comp") or {}).get("chem_comp") or {})
        cid, name = comp.get("id") or "", comp.get("name") or ""
        if not cid:
            continue
        syn = [s.get("name") or "" for s in (n.get("nonpolymer_comp") or {}).get("rcsb_chem_comp_synonyms") or []
               if (s.get("provenance_source") or "") in ("PDB Reference Data", "Author")]
        drug = [s.get("name") or "" for s in (n.get("nonpolymer_comp") or {}).get("rcsb_chem_comp_synonyms") or []
                if (s.get("provenance_source") or "") == "DrugBank"]
        pdesc = ((n.get("rcsb_nonpolymer_entity") or {}).get("pdbx_description") or "")
        named = [w for w in [*COMMON_ALT.get(cid, []), COMMON.get(cid), pdesc, name, *syn, *drug] if w]
        said = _in_title([v for w in named for v in _variants(re.sub(r"<[^>]+>|&\w+;", "", w)) if _readable(v)] +
                         ([cid] if len(cid) >= 3 and not cid.isdigit() else []), title)
        if cid in IONS:
            said = (re.search(rf"(?<![a-z])(?:{ION_WORDS[cid]})(?![\w+-]*\s*(?:-\s*)?{ION_NOUN})", title, re.I)
                            and cid)
        # A lipid by any of its names (the chemical name may be an IUPAC string), except
        # signalling lipids mpstruc names by a short name (PIP<sub>2</sub>).
        quiet = cid in IONS or cid in ADDITIVES or cid in COFACTORS or LIPIDISH.search(name) or \
            (cid not in COMMON and any(LIPIDISH.search(w) for w in [pdesc, *syn]))
        if quiet and not said:
            continue
        how = "title" if said else "name"
        if cid in IONS:
            label = IONS[cid]
        elif said and ((cid in COMMON_ALT and cid not in ("ALF", "BEF", "MGF")) or said.upper() != cid):
            label = said
        elif cid in COMMON:
            label = COMMON[cid]
        else:
            # A readable common name before a systematic one ('auxin'-style synonyms
            # before '1H-indol-3-ylacetic acid'), then the shortest.
            # (DrugBank synonyms include other languages' names, so only the PDB's own are used.)
            nice = list(dict.fromkeys(_bare(_chem_case(w)) for w in [pdesc, name, *syn] if w and _readable(w)))
            if nice and _systematic(nice[0]):
                nice.sort(key=lambda w: (_systematic(w), len(w)))
            if nice:
                label = _chem_case(nice[0])
            else:
                label, how = cid, "id"
        label = _lower_name(label)
        if label not in (l for _, l, _ in out):
            out.append((cid, label, how))
    # Ligands the title names come first, in the title's order.
    pos = lambda lab: (title.lower().find(_plain(lab).lower()) if _plain(lab).lower() in title.lower() else 999)
    out.sort(key=lambda x: pos(x[1]))
    return out[:4]


def _variants(name: str) -> list[str]:
    """Other ways a title writes a chemical name: without stereo prefixes
    ('(S)-(-)-nicotine' -> 'nicotine') and as the anion ('imidodiphosphoric acid'
    -> 'imidodiphosphate'). The last item is the plainest form."""
    out = [name]
    bare = _bare(name)
    if bare != name:
        out.append(bare)
    for w in list(out):
        anion = re.sub(r"(?:(?<=phosph)or|(?<=sulf)ur)?(ic|ous) acid$", lambda m: ("ate" if m[1].lower() == "ic" else "ite") if m[0].islower() else
                       ("ATE" if m[1].upper() == "IC" else "ITE"), w, flags=re.I)
        if anion != w:
            out.append(anion)
    return out


def _bare(name: str) -> str:
    """'(S)-(-)-nicotine' -> 'nicotine', 'D-malate' -> 'malate'."""
    return re.sub(r"^(?:\((?:[RSEZ+-]|\d+[RSEZ])(?:,\s*\d*[RSEZ])*\)-|(?-i:[DL]-(?=[a-z])))+", "", name, flags=re.I) or name


def _systematic(name: str) -> bool:
    """A systematic (IUPAC-style) chemical name rather than a common one."""
    return bool(re.search(r"\d-|-\d|\(\w*\)-|\[|yl[a-z]*[- ][a-z-]*yl|oxy-|amino-", name, re.I))


ELEMENTS = set("H B C N O F P S K V W I U FE MN MG CA CU ZN CO NI MO CL NA SE BR".split())


def _chem_case(name: str) -> str:
    """mpstruc's case for a PDB chemical name written in capitals: drug names in
    lower case ('PAXILLINE' -> 'paxilline') but locants, stereo descriptors and
    formulas kept ('1H-indol-3-ylacetic acid', '(S)-nicotine', 'Fe2/S2 cluster')."""
    if not (name.isupper() and len(name) > 4):
        return name
    words = []
    for w in re.split(r"(\s+)", name):
        letters = re.findall(r"[A-Z]+", w)
        if letters and re.search(r"\d", w) and all(l in ELEMENTS for l in letters):
            words.append(re.sub(r"[A-Z]+", lambda m: m[0][0] + m[0][1:].lower(), w))  # 'FE2/S2' -> 'Fe2/S2'
            continue
        w = w.lower()
        w = re.sub(r"(?<![a-z])(\d+(?:,\d+)*)h-", lambda m: m[1] + "H-", w)          # 1H-, 9H-
        w = re.sub(r"\(([rsez])\)", lambda m: f"({m[1].upper()})", w)                 # (S)-
        w = re.sub(r"\((\d+)([rsez])\)", lambda m: f"({m[1]}{m[2].upper()})", w)
        w = re.sub(r"(?<![a-z])([nos])(,[nos])*-(?=[a-z\d(])", lambda m: m[0].upper(), w)  # N-benzyl, N,N-
        w = re.sub(r"^([dl])-(?=[a-z])", lambda m: m[1].upper() + "-", w)               # L-glutamate
        words.append(w)
    return "".join(words)


def _plain(s: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", "", s))


def _lower_name(label: str) -> str:
    """mpstruc writes drug and toxin names in lower case ('with bound paxilline':
    599 of 697 'with bound' names in 2023-2026 entries); acronyms stay as they are."""
    return label[0].lower() + label[1:] if re.fullmatch(r"[A-Z][a-z]{3,}(?:[ -][A-Za-z0-9]+)*", label) else label


# ---------------------------------------------------------------- what the title says

_SIDE = r"(?:inward|outward|cytoplasmic|periplasmic|extracellular|intracellular|lumenal|cytosolic)"
_DEGREE = r"(?:(?:partially|fully)[- ])?(?:open|closed|occluded)"
STATE_ADJ = (r"apo|open|closed|active|inactive|inactivated|deactive|desensiti[sz]ed|resting|latent|primed|engaged|"
             r"unlatched|latched|pre-?active|pre-?open|pre-?fusion|"
             r"intermediate|occluded|inhibited|auto-?inhibited|ligand-free|nucleotide-free|substrate-free|"
             rf"{_SIDE}[- ](?:facing(?:[- ]{_DEGREE})?|{_DEGREE})|"
             r"(?:partially|fully)[- ](?:open|closed|occluded)|conductive|non-?conductive|post-?translocation|"
             r"pre-?translocation|"
             r"E1|E2|E1P|E2P|E1-ATP|E2-Pi|E2P?-BeF\w*|E1P-ADP|turnover|"
             r"(?:Ca|Na|K|Zn|Mg|Cu|Ca2\+|Na\+|K\+)-(?:bound|free)|"
             r"monomeric|dimeric|trimeric|tetrameric|pentameric|hexameric|(?:\d+-fold[- ])?a?symmetric(?:al)?")
# '4-fold symmetrical open state': a symmetry word qualifying the state it precedes.
STATE_MOD = r"(?:\d+-fold[- ])?a?symmetric(?:al)?"
FACETS = [
    ("ph", re.compile(r"\bpH\s*[-=:]?\s*(\d+(?:\.\d+)?)", re.I), lambda m: f"pH {m[1]}"),
    ("temp", re.compile(r"\b(?:at\s+)?(-?\d{1,3}(?:\.\d)?)\s*(?:°\s*C|degrees? C|deg C|oC)\b", re.I),
     lambda m: f"{m[1]} °C"),
    ("temp", re.compile(r"\b(?:at\s+)(\d{2,3})\s*K\b"), lambda m: f"{m[1]} K"),
    # 'inward facing partially open 2 conformation' -> 'inward facing partially open state 2'
    ("state", re.compile(rf"(?<![\w-])((?:{STATE_MOD}\s+)?(?:{STATE_ADJ}))"
                         r"(?:\s+(\d)(?=[\s-]+(?:state|conformation|form)\b|\s*[,;)]|\s*$))?"
                         r"(?:[\s-]+(?:state|conformation|form|conformational state)(?:\s+(\d|I{1,3}|IV)(?![\w.]))?)?"
                         r"(?![\w-])",
                         re.I), lambda m: _state_phrase(m)),
    ("cls", re.compile(r"\b((?i:state|class|conformation|conformer|form|model|type))(?:\s+(\d{1,2}[a-z]?|[A-D]|I{1,3}[ab]?|IV)|"
                       r"(\d{1,2}))\b(?![-.]\d)"),
     lambda m: f"{m[1].lower()} {m[2] or m[3]}"),
    ("delay", re.compile(r"(?<![\w.])(\d+(?:\.\d+)?\s?(?:fs|ps|ns|µs|μs|us|ms))\b"), lambda m: m[1]),
    # Time points: '300-seconds post reaction initiation' -> '300 seconds'
    ("delay", re.compile(r"(?<![\w.-])(\d+(?:\.\d+)?(?:-\d+(?:\.\d+)?)?)[- ]?(seconds?|secs?|minutes?|mins?|hours?|hrs?)\b",
                         re.I),
     lambda m: f"{m[1]} " + {"sec": "seconds", "secs": "seconds", "min": "minutes", "mins": "minutes", "hr": "hours",
                              "hrs": "hours"}.get(m[2].lower(), m[2].lower())),
    ("dom", re.compile(r"\b(peripheral arm|membrane arm|membrane domain|matrix arm|soluble domain|(?:LBD|NTD|ATD)-TMD|"
                       r"TMD|LBD|NTD|ATD|ECD|ICD|CTD|(?:N|C)-terminal domain|transmembrane domain|extracellular domain|"
                       r"cytoplasmic domain|ligand[- ]binding domain)\b", re.I),
     lambda m: m[1] if m[1].isupper() or "-T" in m[1] else m[1].lower()),
    ("conc", re.compile(r"\b(\d+(?:\.\d+)?)\s?([mµμnu]M)\s+([A-Za-z][\w+-]*)"),
     lambda m: f"{m[1]} {m[2].replace('u', 'µ')} {m[3]}"),
    ("env", re.compile(r"\b(?:in\s+)?(?:(?:lipid|MSP\w*|saposin)\s+)?(nanodiscs?|amphipols?|(?:proteo)?liposomes?|"
                       r"peptidiscs?|SMALPs?|digitonin|LMNG|DDM|GDN|UDM|detergent|bicelles?|"
                       r"native membranes?|membrane vesicles?)\b", re.I),
     lambda m: _env_phrase(m[1])),
    ("cls", re.compile(r"\b(consensus(?: map| refinement)?|composite(?: map| model)?|focused(?: refinement| map)?|"
                       r"local(?:ly)?[- ]refine(?:d|ment)(?: map)?|focus(?:ed)?[- ]refined|(?:gating[- ]ring )?masked|"
                       r"monomer|dimer|trimer|tetramer|pentamer|hexamer|heptamer|octamer|protomer)\b", re.I),
     lambda m: re.sub(r"^gating ring", "gating-ring",
                      re.sub(r"^locally[- ]refined.*|^local[- ]refine.*", "local refinement", m[1].lower()))),
]
# Point mutations: one-letter residue, 2-4 digit position, new residue ('E123A'; 'A2A' is a name, not a mutation).
MUTATION = re.compile(r"(?<![A-Za-z0-9])([ACDEFGHIKLMNPQRSTVWY]\d{2,4}[ACDEFGHIKLMNPQRSTVWY*])(?![A-Za-z0-9])")


def _state_phrase(m: re.Match) -> str:
    adj, num = m[1], (m[2] or m[3]) if m.lastindex and m.lastindex >= 2 else None
    num = f" {num}" if num else ""
    if re.fullmatch(r"apo", adj, re.I):
        return "apo"
    if re.fullmatch(r"(?:mono|di|tri|tetra|penta|hexa)meric", adj, re.I):
        return f"{adj.lower()} form"
    if re.fullmatch(STATE_MOD, adj, re.I):  # '(holo, asymmetric)': curators write just 'asymmetric'
        return adj.lower()
    if re.fullmatch(r"E\w.*", adj) or re.match(r"(?:Ca|Na|K|Zn|Mg|Cu)\b", adj):
        return f"{adj} state{num}"  # 'E2P state', 'Ca2+-bound state'
    return f"{adj.lower()} state{num}"


def _env_phrase(w: str) -> str:
    w = w if w.isupper() else w.lower()
    if re.match(r"nanodisc|amphipol|liposome|proteoliposome|peptidisc|bicelle|smalp", w, re.I):
        return "in " + (w if w.endswith("s") else w + "s")
    return f"in {w}"


def title_facets(title: str) -> dict[str, list[str]]:
    """Condition phrases a PDB title states: pH, temperature, state,
    class/conformation number, membrane mimetic, map type and mutations."""
    out: dict[str, list[str]] = {}
    for kind, rx, fmt in FACETS:
        for m in rx.finditer(title or ""):
            p = fmt(m)
            if kind == "state" and p.split()[0] in ("state", "form"):
                continue
            lst = out.setdefault(kind, [])
            if p not in lst and not any(p in q for q in lst):
                lst.append(p)
    # 'partially open state 1' already says 'state 1'
    if out.get("cls") and out.get("state"):
        out["cls"] = [c for c in out["cls"] if not any(c in s for s in out["state"])]
        if not out["cls"]:
            del out["cls"]
    muts = _dedupe_mutations(MUTATION.findall(title or ""))
    if muts:
        out["mut"] = muts
    return out


def _dedupe_mutations(muts: list[str]) -> list[str]:
    """One substitution per position: 'I4657M' and 'C4657M' cannot both be right,
    so the first one given (the title's, then the depositors') is kept."""
    out, seen = [], set()
    for m in muts:
        pos = m[1:-1]
        if pos not in seen:
            seen.add(pos)
            out.append(m)
    return out


def mutations(det: dict, core: set[str]) -> list[str]:
    """Point mutations the depositors list on the membrane protein's chains."""
    out = []
    for e in _protein_entities(det):
        if core and not (_uniprots(e) & core):
            continue
        out += MUTATION.findall((e.get("rcsb_polymer_entity") or {}).get("pdbx_mutation") or "")
    return _dedupe_mutations(out)


# ---------------------------------------------------------------- one structure's facets

def _aliases(det: dict, core: set[str], extra: list[str] = ()) -> set[str]:
    """Names of the structure's own protein: entity names, gene symbols and the mpstruc name."""
    mine = [e for e in _protein_entities(det) if _uniprots(e) & core] or \
           [e for e in _protein_entities(det) if not _is_antibody(_desc(e))][:1]
    out = set(extra)
    for x in extra:
        out |= set(re.findall(r"\(([^(),]{2,12})\)", x))
    for e in mine:
        out |= set(_genes(e))
        out.add(_desc(e))
        out |= {w for w in re.findall(r"\(([A-Za-z0-9-]{3,})\)", _desc(e))}
    return {a for a in out if a}


# Title words that describe the sample, never a small molecule ('dimer-GPR3-Gs complex').
NOT_A_LIGAND = re.compile(r"^(?:monomers?|dimers?|trimers?|tetramers?|pentamers?|hexamers?|oligomers?|protomers?|"
                          r"heterodimers?|homodimers?|complex(?:es)?|states?|forms?|sugars?|ligands?|substrates?|"
                          r"compounds?|drugs?|inhibitors?|agonists?|antagonists?|lipids?|nascent|apo|holo)$", re.I)


def _chem_names(det: dict, cid: str) -> list[str]:
    """Every name the PDB gives a chemical: its name, description and synonyms."""
    for n in det.get("nonpolymer_entities") or []:
        nc = n.get("nonpolymer_comp") or {}
        if ((nc.get("chem_comp") or {}).get("id") or "") == cid:
            return [w for w in [(nc.get("chem_comp") or {}).get("name"),
                                (n.get("rcsb_nonpolymer_entity") or {}).get("pdbx_description"),
                                *[s.get("name") for s in nc.get("rcsb_chem_comp_synonyms") or []]] if w]
    return []


def _fits_ligand(term: str, names: list[str], own_words: set[str], several: bool) -> bool:
    """Can a title's name `term` be the ligand the PDB knows only by a systematic
    name? It must be one word that is not the protein's own name ('KG-OXGR1') or
    a word for the sample ('dimer'), and look like a short name for a chemical:
    a code or abbreviation ('2-APB', 'UK5099', 'cVA') or a common name that shares
    two stems with one of the PDB's names ('propylamine' for '3-aminopropane').
    With several ligands, an abbreviation's letters must also come in order in
    one of their names ('PHQ' for '2-phenyl-hydroquinone'), so it goes to the right one."""
    if not re.fullmatch(r"\S+", term) or NOT_A_LIGAND.match(term):
        return False
    words = _norm(term).split()
    if any(w in own_words for w in words) or not words:
        return False
    low = [re.sub(r"[^a-z]", "", _plain(n).lower()) for n in names]
    letters = re.sub(r"[^a-z]", "", term.lower())
    if re.search(r"\d", term) and re.fullmatch(r"[A-Za-z]{0,6}-?\d[\w-]*", term) and \
            sum(c.isalpha() for c in term) <= 6:
        return True  # a code name ('UK5099', 'VUF6884', 'G6P', '2-OG', 'F4', '5c')
    if len(re.findall(r"[A-Z]", term)) >= 2 and len(term) <= 8:
        if not several:
            return True
        return any(re.search(".*".join(letters), n) for n in low) if letters else False
    # A common name: two separate 4-letter stems found in one of the PDB's names.
    for n in low:
        hits, i = 0, 0
        while i <= len(letters) - 4:
            if letters[i:i + 4] in n:
                hits, i = hits + 1, i + 4
            else:
                i += 1
        if hits >= 2:
            return True
    return False


def facets(det: dict, core: set[str], use_entity_mutations: bool = True, aliases: list[str] = ()) -> dict:
    """Everything that can tell one structure of a study from another, as
    mpstruc phrases grouped by kind: mut, g (G protein/arrestin), cplx,
    bound, apo, env, ph, temp, state, cls. 'named' lists the partners and
    ligands the PDB title itself names."""
    title = _title(det)
    f: dict = title_facets(title)
    names = set(_desc(e) for e in _protein_entities(det))
    f["mut"] = [m for m in f.get("mut", []) if not any(m in n for n in names | set(aliases))]
    if use_entity_mutations:
        f["mut"] = _dedupe_mutations(f["mut"] + mutations(det, core))
    own = "|".join(_aliases(det, core, aliases))
    pt = partners(det, core, title, own)
    tnames = title_names(title, _aliases(det, core, aliases))
    used = {_norm(x) for x, said in pt["partners"] + pt["peptides"] if said}
    free = [t for t in tnames if _norm(t) not in used]
    ligs = ligands(det, title)
    used |= {_norm(l) for _, l, _ in ligs}
    used |= {_norm(x) for x, _ in pt["partners"]}
    free = [t for t in free if _norm(t) not in used and not any(_norm(t) in u or u in _norm(t) for u in used if u)
            and not any(_norm(t).split()[0] in u.split() for u in used if u and len(_norm(t).split()[0]) >= 4)]
    peptides = []
    # Not the subunits of a large complex (photosystem II, complex I): short chains there are
    # subunits, named only when the title names them.
    n_core = sum(1 for e in _protein_entities(det) if _uniprots(e) & core)
    small = len(pt["partners"]) + len(pt["peptides"]) <= 3 and n_core <= 3
    own_words = set(_norm(own.replace("|", " ")).split())
    for lab, said in pt["peptides"]:
        if not small and not said:
            continue
        if not said and free:
            lab, said = _lower_name(free.pop(0)), True
        peptides.append((lab, said))
    bound = [(lab, how) for _, lab, how in ligs]
    # A small molecule known only by its chemical ID, then one known only by a
    # systematic name ('2-aminoethyl diphenylborinate'), takes the name the title
    # gives it ('2-APB'; the systematic case is from PR #4), unless that name is a
    # word of the protein's own name ('A2A').
    order = [i for i, b in enumerate(bound) if b[1] == "id"] + \
            [i for i, b in enumerate(bound) if b[1] == "name" and _systematic(b[0])]
    for i in order:
        if bound[i][1] == "id":
            pick = next((t for t in free if _norm(t) not in own_words and not NOT_A_LIGAND.match(t) and
                         re.fullmatch(r"\S+(?: [a-z]\S*)?(?: [A-Z0-9]{1,3})?", t)), None)
        else:
            names = _chem_names(det, ligs[i][0]) + [bound[i][0]]
            pick = next((t for t in free if _fits_ligand(t, names, own_words, len(ligs) > 1)), None)
        if pick:
            free.remove(pick)
            bound[i] = (_lower_name(pick), "title")
    bound = [b for b in bound if b[1] != "id"]  # a bare chemical ID ('A1EYI') means nothing to a reader
    # One ligand the title does not name and one compound name the title gives: the
    # title's name is the one curators use ('LSD', not 'lysergic acid diethylamide').
    unnamed = [i for i, b in enumerate(bound) if b[1] == "name"]
    spare = [t for t in free if _norm(t) not in own_words and re.fullmatch(r"\S+(?: [A-Z0-9])?", t)]
    if len(unnamed) == 1 and len(spare) == 1 and len(bound[unnamed[0]][0].split()) >= 3:
        bound[unnamed[0]] = (_lower_name(spare[0]), "title")
    bound = [(lab, how == "title") for lab, how in bound]
    f["g"] = [pt["g"]] if pt["g"] else []
    f["g_from"] = pt["g_from"]
    f["cplx"] = [x for x, _ in pt["partners"]][:3] if small else [x for x, s in pt["partners"] if s][:3]
    f["bound"] = [x for x, _ in peptides] + [x for x, _ in bound]
    f["named"] = {x for x, s in pt["partners"] + peptides + bound if s}
    f["ab"] = [x for x, s in pt["antibodies"] if s]
    if "apo" in f.get("state", []):
        f["state"].remove("apo")
        f["apo"] = ["apo"]
    if f.get("conc"):
        said_conc = {_norm(c.split()[-1]) for c in f["conc"]}
        f["bound"] = [b for b in f["bound"] if _norm(b) not in said_conc]
    # 'in GDN detergent': the detergent's name is enough
    if len(f.get("env", [])) > 1 and "in detergent" in f["env"]:
        f["env"].remove("in detergent")
    return f


def _states(states: list[str]) -> list[str]:
    """'inward-open state', 'substrate-free state' -> 'inward-open, substrate-free state'."""
    if len(states) < 2 or not all(x.endswith(" state") for x in states):
        return list(states)
    return [", ".join(x[:-len(" state")] for x in states) + " state"]


def _mut_phrase(muts: list[str]) -> str:
    return "/".join(muts) + " mutant"


def _phrase(f: dict, method_prefix: str = "") -> str:
    """Facets -> one mpstruc phrase, in curator order."""
    bits = []
    if f.get("mut"):
        bits.append(_mut_phrase(f["mut"]))
    lead = ""
    if f.get("g"):
        lead = f"{f['g'][0]} complex"
    if f.get("cplx"):
        lead = (lead + " " if lead else "") + "in complex with " + _join(f["cplx"])
    if f.get("bound"):
        lead = (lead + " " if lead else "") + "with bound " + _join(f["bound"])
    elif f.get("apo"):
        lead = (lead + ", " if lead else "") + "apo"
    if f.get("conc"):
        lead = (lead + ", " if lead else "") + "with " + _join(f["conc"])
    if lead:
        bits.append(lead)
    if f.get("dom"):
        bits.insert(0, ", ".join(f["dom"]))
    for k in ("env", "ph", "temp", "delay", "state", "cls", "ab", "rest"):
        bits += [f"in complex with {x}" for x in f.get(k) or []] if k == "ab" else \
            _states(f.get(k) or []) if k == "state" else f.get(k) or []
    s = ", ".join(bits)
    return f"{method_prefix}{s}" if s else method_prefix.rstrip(": ")


# ---------------------------------------------------------------- what differs within a study

def _title_words(title: str) -> list[str]:
    # 'in situ' stays one word; '2:2' keeps its colon.
    t = re.sub(r"\b(in|ex) (situ|vivo|vitro|silico)\b", r"\1~\2", title or "", flags=re.I)
    return [w.replace("~", " ") for w in re.findall(r"[A-Za-z0-9α-ωΑ-Ω][\w+().'/~-]*(?::\d+)*", t)]


TITLE_NOISE = re.compile(r"^(?:cryo-?em|cryoem|structure|structures|of|human|the|a|an|determined|obtained|by|"
                         r"cryo-electron|microscopy|x-ray|crystal|resolution|angstrom|angstroms|å|a|map|model|"
                         r"reconstruction|data|dataset|collected|using)$", re.I)
# Words that say nothing on their own: dropped from the ends of a run of title words.
FILLER = re.compile(r"^(?:presence|absence|complex|complexed|state|structure|bound|form|conformation|"
                    r"basal|condition|conditions|sample|post|after|before|following)$", re.I)
NUMBER_UNIT = re.compile(r"^(?:kev|ev|å|a|angstroms?|kda|da|mm|µm|μm|nm|k)$", re.I)


def _distinct_title(titles: list[str], i: int, covered: str) -> str | None:
    """Fallback: the run(s) of words in title i that the other titles lack.
    Decimal numbers (resolutions, '9.5 keV'), numbers with units and filler words
    are left out; whole numbers that name a class or state ('visible 2') stay."""
    key = lambda w: w.lower().strip("().,")
    sets = [{key(w) for w in _title_words(t)} for t in titles]
    common = set.intersection(*sets) if sets else set()
    words = _title_words(titles[i])
    runs, run = [], []
    taken = -1
    for j, w in enumerate(words):
        if j == taken:
            continue
        k = key(w)
        nxt = key(words[j + 1]) if j + 1 < len(words) else ""
        if re.fullmatch(r"\d+\.\d+", k) or (re.fullmatch(r"\d+", k) and NUMBER_UNIT.match(nxt)) or \
                (NUMBER_UNIT.match(k) and j and re.fullmatch(r"\d+(?:\.\d+)?", key(words[j - 1]))):
            if run:  # a measurement ends the run and is dropped, with its unit
                runs.append(run)
            run = []
            continue
        ident = j > 0 and re.fullmatch(r"[A-Z]", w)  # 'chain A', 'form B'
        if k not in common and (ident or (k not in STOP and not TITLE_NOISE.match(k))):
            run.append(w.strip("(),"))
            if re.fullmatch(r"\d+", k) and len(run) == 1 and nxt and nxt not in STOP:
                run.append(words[j + 1].strip("(),"))  # '2493 GW/cm2', not a bare '2493'
                taken = j + 1
        elif run and k in STOP and k not in common:
            run.append(w)
        else:
            if run:
                runs.append(run)
            run = []
    if run:
        runs.append(run)
    cov = _norm(covered)
    out = []
    for r in runs:
        while r and (r[-1].lower() in STOP or (FILLER.match(r[-1]) and len(r) == 1)):
            r = r[:-1]
        while r and (r[0].lower() in STOP or FILLER.match(r[0])):
            r = r[1:]
        txt = " ".join(r)
        if txt and _norm(txt) not in cov and not re.fullmatch(r"[\d.]+", txt) and txt not in out:
            out.append(txt)
    s = ", ".join(out)
    return s if 0 < len(s.split()) <= 8 else None


PRIMARY = ("mut", "g", "cplx", "bound", "conc", "ph", "state", "dom", "ab")  # compared with the whole group
MINOR = ("temp", "env", "cls", "delay")  # shown only when the related structures differ among themselves


def differences(group: list[tuple[str, dict | None]], core: set[str], aliases: list[str] = ()) -> dict[str, str]:
    """{pdb: phrase} naming only what sets each structure of `group` apart
    (the first structure is the lead, described by the entry's name).

    Whatever every structure shares (bound partners, ligands, pH, mutations,
    membrane mimetic ...) is left out; if only the ligands differ the phrase
    is just 'with bound X'. When nothing in the PDB data differs, the words
    of its title that the other titles lack are used, then 'additional
    structure'."""
    present = [(p, d) for p, d in group if d]
    if not present:
        return {p: "additional structure" for p, _ in group}
    titles = [_title(d) for _, d in present]
    # Engineering mutations listed on every chain are noise unless a title speaks of mutants.
    talk_mut = any(re.search(r"mutant|mutation|variant|substitution", t, re.I) or MUTATION.search(t) for t in titles)
    fs = [facets(d, core, use_entity_mutations=talk_mut, aliases=aliases) for _, d in present]
    gkey = _g_keys(fs)
    methods = [_method([m.get("method", "") for m in d.get("exptl") or []]) for _, d in present]
    others = fs[1:] if len(fs) > 2 else fs
    out = {}
    for i, (pdb, d) in enumerate(present):
        own: dict = {}
        for k in PRIMARY + MINOR:
            mine = fs[i].get(k) or []
            pool = fs if k in PRIMARY or i == 0 else others
            kf = gkey if k == "g" else (lambda x: x)
            common = set.intersection(*({kf(x) for x in f.get(k) or []} for f in pool))
            if k in PRIMARY and i > 0 and len(fs) > 3:
                # What the lead has and most of the study shares goes without saying
                # ('latent state', not 'G<sub>i</sub> complex, latent state', when one
                # structure of the study is a nanobody complex instead).
                lead = {kf(x) for x in fs[0].get(k) or []}
                common |= {x for x in lead if sum(x in {kf(y) for y in f.get(k) or []} for f in fs) * 2 > len(fs)}
            keep = [x for x in mine if kf(x) not in common]
            if keep:
                own[k] = keep
        if fs[i].get("apo") and not all(f.get("apo") for f in fs):
            own["apo"] = ["apo"]
        prefix = ""
        if methods[i] and methods[0] and methods[i] != methods[0]:
            prefix = methods[i].replace(" structure", "") + ": "
        phrase = _phrase(own, prefix)
        if not phrase or phrase == prefix.rstrip(": "):
            rest = _distinct_title(titles, i, phrase)
            if not rest and i > 0 and set(fs[i].get("bound") or []) != set(fs[0].get("bound") or []):
                # It only lacks a ligand the lead has: say what it does have.
                rest = _phrase({"bound": fs[i]["bound"]} if fs[i].get("bound") else {"apo": ["apo"]})
            phrase = f"{prefix}{rest}" if rest else (prefix + "additional structure")
        out[pdb] = phrase
    _tell_apart(out, present, titles)
    for p, _ in group:
        out.setdefault(p, "additional structure")
    return out


def _g_keys(fs: list[dict]):
    """How G proteins are compared within a study: 'G<sub>i1</sub>' is the same
    partner as 'G<sub>i</sub>'; and when the titles name one subtype ('Gq') but
    the chain names give another ('Gi', the same engineered chimera), they are
    one partner too, so neither is listed as a difference."""
    base = lambda x: re.sub(r"(<sub>[a-z]+)\d+(</sub>)", r"\1\2", x)
    alias: dict[str, str] = {}
    gs = [(f["g"][0], f.get("g_from")) for f in fs if f.get("g")]
    if len(gs) == len(fs) and not any("arrestin" in g for g, _ in gs):
        by_title = {base(g) for g, src in gs if src == "title"}
        by_chain = {base(g) for g, src in gs if src != "title"}
        if len(by_title) == 1 and len(by_chain) == 1 and by_title != by_chain:
            alias[next(iter(by_chain))] = next(iter(by_title))
    return lambda x: alias.get(base(x), base(x))


def _tell_apart(out: dict[str, str], present: list[tuple[str, dict]], titles: list[str]) -> None:
    """Structures of a study whose phrases came out the same although their PDB
    titles differ get the title words that set them apart from each other
    ('open state, 4-fold symmetrical' vs 'open state, 2-fold symmetrical')."""
    groups: dict[str, list[int]] = {}
    for i, (pdb, _) in enumerate(present):
        groups.setdefault(out[pdb], []).append(i)
    for phrase, idx in groups.items():
        tt = [titles[i] for i in idx]
        if len(idx) < 2 or len({_norm(t) for t in tt}) < 2:
            continue
        for n, i in enumerate(idx):
            extra = _distinct_title(tt, n, phrase if phrase != "additional structure" else "")
            if not extra:
                continue
            pdb = present[i][0]
            out[pdb] = phrase[:-len("additional structure")] + extra if phrase.endswith("additional structure") \
                else f"{phrase}, {extra}"


# ---------------------------------------------------------------- names

_COND = re.compile(
    r"(?:,\s+|\s+-\s+|\s*-\s*(?=(?:mini[- ]?)?G<sub>|(?:mini[- ]?)?G\s?protein|β-arrestin|&beta;-arrestin|beta-arrestin|"
    r"(?:mini[- ]?)?G[soqi]\b)|\s+|;\s*)"
    r"(?=with\b|bound\b|in complex\b|complexed\b|apo\b|in (?:the|an?)\b|in (?-i:[A-Z][A-Z0-9]{2,})\b|in (?:lipid|nanodisc|detergent|LMNG|DDM|GDN|"
    r"digitonin|amphipol|MSP|saposin|peptidisc|total|native|liposome|presence|absence)|in [a-z]+[- ]facing|"
    r"[A-Z]\d+[A-Z](?:/[A-Z]\d+[A-Z])* mutant|mutant\b|[\w-]+-(?:bound|free)\b|[\w-]+[- ]facing\b|"
    r"(?:open|closed|inactive|active|resting|desensiti[sz]ed|occluded|inhibited|intermediate)\b|pH\b|E\d\b|"
    r"state\b|form\b|(?:mini[- ]?)?G<sub>|G protein|β-arrestin|&beta;-arrestin|beta-arrestin|focused|local\b|"
    r"consensus|composite|full[- ]length|crystal|x-ray|cryo|structure\b)", re.I)


def base_name(name: str | None) -> str | None:
    """An mpstruc name without its condition:
    'C5a receptor (C5aR) - G<sub>i</sub> complex, with bound C5a' -> 'C5a receptor (C5aR)'.
    Commas inside brackets ('(MCT2, SLC16 family)') or words ('Inositol-1,4,5-')
    do not end the protein name."""
    if not name:
        return None
    s = re.sub(r"\s+", " ", name).strip()
    cut = len(s)
    for m in re.finditer(r",\s|;\s*", s):
        # An entity's ';' ('&alpha;', '&Delta;') does not end the name.
        if m[0][0] == ";" and re.search(r"&#?\w+$", s[:m.start()]):
            continue
        if _depth_at(s, m.start()) == 0:
            cut = m.start()
            break
    for m in _COND.finditer(s):
        if m.start() < cut and _depth_at(s, m.start()) == 0 and m.start() > 0:
            # '&quot; with ...': the ';' belongs to the entity.
            cut = m.start() + (s[m.start()] == ";" and bool(re.search(r"&#?\w+$", s[:m.start()])))
            break
    out = s[:cut]
    # Lower case only: 'Monoamine Oxidase A' keeps its 'A'. A word that only led into
    # the setting goes too ('RyR1 ryanodine receptor embedded in a lipid bilayer').
    out = re.sub(r"(?:\s+(?:in|with|and|of|at|from|by|for|to|the|a|an|embedded|reconstituted|incorporated|"
                 r"solubili[sz]ed|inserted|expressed|purified|captured|trapped|determined|bound))+$", "", out)
    out = out.strip(" ,-") if re.search(r"&#?\w+;$", out.strip(" ,-")) else out.strip(" ,;-")
    return out or name.strip()


def _depth_at(s: str, i: int) -> int:
    d = 0
    for ch in s[:i]:
        d += ch in "(["
        d -= ch in ")]"
    return max(d, 0)


def family_base(names: list[tuple[str | None, str | None]]) -> str | None:
    """The protein name mpstruc uses most for a family, from (name, year)
    pairs of the entries that share the protein; ties go to the newest.
    Names that differ only in case or markup ('D3' / 'D<sub>3</sub>') count
    as one name, written in its most common spelling (PR #4 listed picking
    the newest spelling as a known gap)."""
    seen: dict[str, dict] = {}
    for name, year in names:
        b = base_name(name)
        if not b:
            continue
        rec = seen.setdefault(_norm(b), {"n": 0, "year": "", "spell": Counter(), "last": {}})
        rec["n"] += 1
        rec["year"] = max(rec["year"], year or "")
        rec["spell"][b] += 1
        rec["last"][b] = max(rec["last"].get(b, ""), year or "")
    if not seen:
        return None
    rec = max(seen.values(), key=lambda r: (r["n"], r["year"]))
    best = max(rec["spell"], key=lambda b: (rec["spell"][b], rec["last"][b]))
    return _cap(best) if re.match(r"[a-z]{4,}\b", best) else best


# Words that do not make a protein name on their own.
_GENERIC = set("""protein proteins subunit subunits chain chains polypeptide domain family putative
uncharacterized transmembrane membrane complex fragment precursor probable isoform like""".split())
_EDGE = set("""of from the a an in with and bound by to at for on structure""".split())


def complex_name(title: str | None, mine: list[dict], known: list[str] = ()) -> str | None:
    """The name of a complex of several proteins, as its PDB title writes it: the
    longest title phrase that enough of its chains' names contain ('Toxoplasma gondii
    cytochrome c oxidase from the respiratory supercomplex III2-IV', with chains
    'Cytochrome c oxidase subunit 3', 'Cytochrome c oxidase polypeptide II', ... ->
    'cytochrome c oxidase'). None unless at least two chains and a third of them
    share it. When an mpstruc name in `known` (the subfamily's parents) is the same
    phrase, its spelling is used ('Cytochrome C Oxidase')."""
    names = list(dict.fromkeys(_desc(e).lower() for e in mine if _desc(e)))
    words = re.findall(r"[A-Za-z0-9][\w()'+-]*", title or "")
    best = None
    for n in range(min(6, len(words)), 1, -1):  # one word ('Flagellar', 'ATP') is not a name
        for i in range(len(words) - n + 1):
            ws = words[i:i + n]
            if ws[0].lower() in _EDGE or ws[-1].lower() in _EDGE:
                continue
            if len(ws[-1]) == 1 and ws[-1] != "I":  # 'cytochrome c' of 'cytochrome c reductase'
                continue
            if not any(w.lower() not in _GENERIC and w.lower() not in _EDGE and len(w) > 2 for w in ws):
                continue
            phrase = " ".join(ws)
            pat = re.compile(r"(?<![\w-])" + re.escape(phrase.lower()) + r"(?![\w-])")
            hits = sum(bool(pat.search(x)) for x in names)
            if hits >= 2 and 3 * hits >= len(names) and (best is None or (n, hits) > (best[1], best[0])):
                best = (hits, n, phrase)
    if not best:
        return None
    phrase = best[2]
    for k in known:
        b = base_name(k)
        if b and _norm(re.sub(r"<[^>]+>", "", b)) == _norm(phrase):
            return b
    return _cap(phrase)


def pdb_name(mine: list[dict]) -> str | None:
    """Protein name from the PDB entity names, mpstruc-style: 'SLC16A7
    Monocarboxylate transporter 2' (gene symbol first when the name does not
    already contain it, as mpstruc does for new proteins)."""
    names = []
    for e in mine:
        if SEQUENCE.fullmatch(_desc(e).upper()):
            continue  # a bare residue sequence is not a protein name
        d = _clean_name(_desc(e))
        # 'Calcium-transporting ATPase 8, plasma membrane-type' -> 'Calcium-transporting ATPase 8'
        d = re.sub(r",\s+(?=[a-z])[^,]*$", "", d)
        if d and d not in names:
            names.append(d)
    if not names:
        return None
    name = " / ".join(names[:2])
    # PDB charges as mpstruc writes them: 'H(+)/Cl(-) exchange' -> 'H<sup>+</sup>/Cl<sup>-</sup> exchange'
    name = re.sub(r"\b([A-Z][a-z]?)\((\d?[+-])\)", r"\1<sup>\2</sup>", name)
    # The gene symbol only goes before a single protein name, never before a list of fused parts.
    # When the PDB lists a fusion tag's gene first (cybC for BRIL, GFP), no gene is added, as
    # before the tag was dropped from the name: the protein's own gene would then go before
    # receptor names the PDB already spells out ('HTR1A 5-hydroxytryptamine receptor 1A'),
    # which curators do not write, and it scored no better against their names.
    genes = _genes(mine[0]) if len(names) == 1 and ", " not in name else []
    gene = genes[0] if genes and genes[0].lower() not in FUSION_TAG_GENES else None
    if gene and re.fullmatch(r"[A-Za-z][A-Za-z0-9-]{1,9}", gene) and re.search(r"\d|^[A-Z]{3,}$", gene) \
            and not re.fullmatch(r"[A-Z]{1,4}\d{3,}\w*|[a-z]{2,}\d{3,}\w*", gene) \
            and not any(_norm(g) in _norm(name).split() for g in genes):
        name = f"{gene} {name}"
    return _cap(name)


def name_condition(det: dict, core: set[str], aliases: list[str] = ()) -> str:
    """The lead structure's own condition as it goes into its mpstruc name:
    ' - G<sub>s</sub> complex with bound X', ' in complex with Y', ', apo form'.
    Partners and ligands are named only when the PDB title names them or
    they are drugs/substrates; complex subunits and cofactors are not."""
    f = facets(det, core, use_entity_mutations=False, aliases=aliases)
    named = f.get("named", set())
    bits = ""
    if f.get("mut"):
        bits += ", " + _mut_phrase(f["mut"])
    if f.get("g"):
        bits += f" - {f['g'][0]} complex"
    cplx = [x for x in f.get("cplx", []) if x in named]
    if cplx:
        bits += " in complex with " + _join(cplx)
    if f.get("bound"):
        bits += " with bound " + _join(f["bound"])
    tail = []
    if f.get("apo") and not f.get("bound"):
        tail.append("apo form")
    tail += f.get("env", []) + f.get("ph", []) + _states(f.get("state", []))
    if tail:
        bits += ", " + ", ".join(tail)
    return bits


# ---------------------------------------------------------------- the draft

def build(det: dict | None, core: set[str], related: list[tuple[str, dict | None]] = (),
          context: dict | None = None) -> dict | None:
    """Draft mpstruc fields for one structure.

    `related` lists (pdb, details) of structures from the same study that are
    suggested as its related entries; each becomes a description line naming
    only what differs ('pH 6.0, deactive state, 3.10 &Aring;: 9R13').

    `context` (optional) carries what mpstruc already holds:
      family_names  [(name, year)] of mpstruc entries with the same protein:
                    their protein name is reused instead of the PDB's;
      subgroup_names  names of the parents in the suggested subfamily: a
                    complex named from its title takes their spelling;
      group         [(pdb, details)] of the entry this one would be listed
                    under and its related structures: used to write
                    `relatedLine`, the line to add to that entry.
    Text fields are mpstruc raw markup (entities escaped; <br />, <sub>,
    <sup> allowed).
    """
    if not det:
        return None
    context = context or {}
    ents = _protein_entities(det)
    mine = [e for e in ents if _uniprots(e) & core]
    mine = mine or [e for e in ents if not _is_antibody(_desc(e))] or ents
    orgs = [o for e in mine for o in e.get("rcsb_entity_source_organism") or []]
    hosts = [o.get("ncbi_scientific_name") for e in mine for o in e.get("rcsb_entity_host_organism") or []]
    species = next((o.get("ncbi_scientific_name") for o in orgs if o.get("ncbi_scientific_name")), None)
    lineage = {l.get("name") for o in orgs for l in o.get("taxonomy_lineage") or []}
    domain = next((d for d in DOMAINS if d in lineage), None)
    res = _res(det)

    fam = family_base(context.get("family_names") or [])
    base = fam or pdb_name(mine)
    # A complex of several proteins (or chains with only generic names) goes by the
    # complex's name in the title, not its first two chains ('Transmembrane protein /
    # Cg8 family protein' for 9I4Y, a cytochrome c oxidase).
    named = list(dict.fromkeys(_desc(e) for e in mine if _desc(e)))
    plain = base and not fam and not (set(_norm(base).split()) - _GENERIC)
    cplx = None
    if not fam and (len(named) >= 3 or plain):
        cplx = complex_name((det.get("struct") or {}).get("title"), mine, context.get("subgroup_names") or [])
    if cplx:
        base, fam = cplx, cplx if "<" in cplx or "&" in cplx else None
    aliases = [html.unescape(re.sub(r"<[^>]+>", "", base))] if base else []
    name = None
    if base:
        # A reused mpstruc name is already raw markup; a PDB name is plain text.
        name = (base if fam else _mp(base)) + _mp(name_condition(det, core, aliases))
        name = format_text(name)

    lines = []
    method = _method([m.get("method", "") for m in det.get("exptl") or []])
    if method:
        lines.append(method)
    pdb = (det.get("rcsb_id") or "").upper() or "_LEAD"
    if related:
        diff = differences([(pdb, det)] + list(related), core, aliases)
        for rp, rd in related:
            rres = _res(rd)
            lines.append(f"{_mp(diff[rp])}{f', {rres:.2f} &Aring;' if rres else ''}: {rp}")
    related_line = None
    group = [g for g in context.get("group") or [] if g[0] != pdb]
    if group:
        diff = differences(group + [(pdb, det)], core, aliases)
        related_line = f"{_mp(diff[pdb])}{f', {res:.2f} &Aring;' if res else ''}: {pdb}"

    cit = det.get("rcsb_primary_citation") or {}
    pm = cit.get("pdbx_database_id_PubMed")
    pages = "-".join(dict.fromkeys(p for p in (cit.get("page_first"), cit.get("page_last")) if p)) or None
    out = {
        "name": name,
        "species": _esc(species) if species else None,
        "taxonomicDomain": domain,
        "expressedInSpecies": _esc(hosts[0]) if hosts and hosts[0] else None,
        "resolution": f"{res:.2f}".rstrip("0").rstrip(".") if res else None,
        "description": "<br /> ".join(lines) or None,
        "bibliography": {
            "pubMedId": str(pm) if pm and int(pm) > 0 else None,
            "authors": html.escape(_authors(cit.get("rcsb_authors")) or "", quote=False) or None,
            "year": str(cit["year"]) if cit.get("year") else None,
            # As the journal printed it: mpstruc adds no markup to 2024-2026 paper titles.
            "title": html.escape(cit["title"], quote=False) if cit.get("title") else None,
            "journal": html.escape(cit["journal_abbrev"], quote=False)
            if cit.get("journal_abbrev") and cit["journal_abbrev"] != "To Be Published" else None,
            "volume": cit.get("journal_volume"), "issue": cit.get("journal_issue"), "pages": pages,
            "doi": cit.get("pdbx_database_id_DOI"),
        },
        "pdbTitle": _title(det) or None,
    }
    if related_line:
        out["relatedLine"] = related_line
    return out


# ---------------------------------------------------------------- extra PDB fields

EXTRA_QUERY = (
    "query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id "
    "polymer_entities { rcsb_polymer_entity { pdbx_description pdbx_mutation } "
    "entity_poly { rcsb_sample_sequence_length } "
    "rcsb_entity_source_organism { rcsb_gene_name { value } } } "
    "nonpolymer_entities { rcsb_nonpolymer_entity { pdbx_description } "
    "nonpolymer_comp { chem_comp { id } rcsb_chem_comp_synonyms { name provenance_source } } } } }")


def enrich(details: dict[str, dict | None], batch: int = 100, max_age_s: float = 30 * 86400) -> dict[str, dict | None]:
    """Add the fields a listing needs beyond rcsb.entry_details (mutations,
    chain lengths, gene names, ligand synonyms) to each entry, in place.
    Cached in cache/listing_extras.json; skipped quietly if the PDB is unreachable."""
    from . import rcsb
    path = rcsb.CACHE_DIR / "listing_extras.json"
    cache: dict = {}
    if path.exists() and time.time() - path.stat().st_mtime < max_age_s:
        cache = json.loads(path.read_text())
    todo = sorted({p.upper() for p, d in details.items() if d} - cache.keys())
    try:
        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            data = rcsb._post_json(rcsb.GRAPHQL, {"query": EXTRA_QUERY, "variables": {"ids": chunk}})
            found = {e["rcsb_id"]: e for e in (data.get("data") or {}).get("entries") or [] if e}
            for pid in chunk:
                cache[pid] = found.get(pid)
        if todo:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    except Exception as err:  # the listing still works without these fields
        print(f"listing: could not fetch extra PDB fields ({err})")
    for pid, det in details.items():
        x = cache.get(pid.upper()) if det else None
        if x:
            details[pid] = det = copy.deepcopy(det)  # never alter rcsb's cached records
            _merge(det, x)
    return details


def _merge(det: dict, x: dict) -> None:
    pool = list(x.get("polymer_entities") or [])
    for e in det.get("polymer_entities") or []:
        d = _desc(e)
        hit = next((p for p in pool if _desc(p) == d), None)
        if not hit:
            continue
        pool.remove(hit)
        e.setdefault("rcsb_polymer_entity", {})["pdbx_mutation"] = (hit.get("rcsb_polymer_entity") or {}).get("pdbx_mutation")
        e.setdefault("entity_poly", {})["rcsb_sample_sequence_length"] = (hit.get("entity_poly") or {}).get("rcsb_sample_sequence_length")
        genes = [o.get("rcsb_gene_name") for o in hit.get("rcsb_entity_source_organism") or []]
        for o, g in zip(e.get("rcsb_entity_source_organism") or [], genes):
            o["rcsb_gene_name"] = g
    syn = {((n.get("nonpolymer_comp") or {}).get("chem_comp") or {}).get("id"): n for n in x.get("nonpolymer_entities") or []}
    for n in det.get("nonpolymer_entities") or []:
        cid = ((n.get("nonpolymer_comp") or {}).get("chem_comp") or {}).get("id")
        if cid in syn:
            n["rcsb_nonpolymer_entity"] = syn[cid].get("rcsb_nonpolymer_entity")
            n["nonpolymer_comp"]["rcsb_chem_comp_synonyms"] = (syn[cid].get("nonpolymer_comp") or {}).get("rcsb_chem_comp_synonyms")
