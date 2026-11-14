"""Which bound ligands do mpstruc curators name? Learned from mpstruc itself.

A PDB entry lists every non-polymer molecule in the model: the drug the paper
is about, but also detergent, lipid fragments, buffer and ions. Curators name
only some of them ('with bound ATP'), and which ones is not a fixed list. This
module measures it:

* LEARN (`python -m mpstruc.ligands build`). For every structure mpstruc lists
  (entries and related structures) it takes the PDB's non-polymer ligands and
  checks whether the curators' own text for that structure (the condition
  part of the entry name, its description, the related line) mentions each
  one. Per chemical component and per ligand class this gives 'seen in N
  structures, named in M', split by whether the PDB title mentions the ligand.
  The result is mpstruc/ligand_usage.json.
* DECIDE (`decide`). For a ligand of a new structure: the rate at which
  curators named this component (or, for a component mpstruc has rarely or
  never seen, its class) in the same situation (title mentions it / does
  not), compared with a threshold chosen on tuning families. Every decision
  comes with a plain-English reason for the curators. 'The title mentions
  it' is listing.ligands_in_title, the same test in learning and deciding.
  (The paper's title was tried as a second place to look; it did not help.)
* MEASURE (`tune`, `score`). Families of entries are split 60/20/20: usage is
  learned on the first part, the rule chosen on the second, and precision
  and recall of 'ligand named' reported on the third, which nothing saw.

Classes come from the chemical component's own properties (formula, weight,
type, SMILES) and name patterns, never from a list of 'known' ligands, so a
component nobody has seen before still gets a class.
"""
from __future__ import annotations

import html
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
USAGE_PATH = ROOT / "ligand_usage.json"

# Never ligands, whatever the data says: water and placeholders for unknown atoms.
HARD_EXCLUDE = {"HOH", "DOD", "UNX", "UNL", "UNK"}


# ---------------------------------------------------------------- PDB data (cached)

CORPUS_QUERY = (
    "query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id struct { title } "
    "rcsb_primary_citation { title } "
    "polymer_entities { rcsb_polymer_entity { pdbx_description } } "
    "nonpolymer_entities { rcsb_nonpolymer_entity { pdbx_description } "
    "nonpolymer_comp { chem_comp { id name } } } } }")
COMP_QUERY = (
    "query($ids: [String!]!) { chem_comps(comp_ids: $ids) { chem_comp { id name formula formula_weight type } "
    "rcsb_chem_comp_synonyms { name provenance_source } "
    "drugbank { drugbank_info { name } } "
    "rcsb_chem_comp_descriptor { SMILES } } }")


def _cached_fetch(cache_name: str, ids: list[str], query: str, field: str, key, batch: int,
                  slim=lambda x: x) -> dict[str, dict | None]:
    """ids -> records from the PDB's GraphQL API, kept in cache/<cache_name>
    (a component or a released entry never changes in the ways used here)."""
    from . import rcsb
    path = rcsb.CACHE_DIR / cache_name
    cache: dict = json.loads(path.read_text()) if path.exists() else {}
    todo = sorted({i.upper() for i in ids if i} - cache.keys())
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        data = rcsb._post_json(rcsb.GRAPHQL, {"query": query, "variables": {"ids": chunk}})
        found = {key(e): slim(e) for e in (data.get("data") or {}).get(field) or [] if e}
        for cid in chunk:
            cache[cid] = found.get(cid)
        if (i // batch) % 10 == 9 or i + batch >= len(todo):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    return {i: cache.get(i.upper()) for i in ids}


def _slim_comp(e: dict) -> dict:
    c = e.get("chem_comp") or {}
    syn = [[s.get("name"), s.get("provenance_source")] for s in e.get("rcsb_chem_comp_synonyms") or [] if s.get("name")]
    return {"id": c.get("id"), "name": c.get("name"), "formula": c.get("formula"),
            "weight": c.get("formula_weight"), "type": c.get("type"), "synonyms": syn,
            "drug": ((e.get("drugbank") or {}).get("drugbank_info") or {}).get("name"),
            "smiles": (e.get("rcsb_chem_comp_descriptor") or {}).get("SMILES")}


def chem_comps(ids: list[str], batch: int = 150) -> dict[str, dict | None]:
    """Chemical component records {id, name, formula, weight, type, synonyms
    [[name, source]], drug (DrugBank name), smiles}; cached in cache/chem_comps.json."""
    return _cached_fetch("chem_comps.json", ids, COMP_QUERY, "chem_comps",
                         lambda e: (e.get("chem_comp") or {}).get("id"), batch, _slim_comp)


def corpus_entries(ids: list[str], batch: int = 100) -> dict[str, dict | None]:
    """Title, citation title, protein names and non-polymer ligands per PDB entry;
    cached in cache/ligand_corpus.json."""
    return _cached_fetch("ligand_corpus.json", ids, CORPUS_QUERY, "entries", lambda e: e.get("rcsb_id"), batch)


# ---------------------------------------------------------------- ligand classes

CLASS_LABEL = {
    "ion": "single-atom ion", "cofactor": "cofactor, pigment or metal cluster", "nucleotide": "nucleotide",
    "detergent": "detergent", "lipid": "lipid or fatty acyl chain", "sugar": "sugar",
    "additive": "buffer, cryoprotectant or precipitant", "amino acid": "amino acid",
    "drug": "small-molecule ligand"}
CLASSES = tuple(CLASS_LABEL)

_METALS = set("Fe Cu Mo Ni Co Mn Mg Zn W V Cd Hg Pt Au Ag Ru Rh Os Ir Pd Gd Tb Yb Eu Sm La Ca Ba Sr".split())
_ELEMENT_ION = re.compile(
    r"^(?:lithium|sodium|potassium|rubidium|c[ae]?esium|magnesium|calcium|strontium|barium|manganese|iron|fe|cobalt|"
    r"nickel|copper|zinc|cadmium|mercury|chloride|bromide|iodide|fluoride|thallium|gadolinium|ytterbium|platinum|"
    r"gold|silver|lead|xenon|krypton|argon|aluminum|aluminium)\b[ (]*(?:[IV]+|\d?[+-])?\)? ?(?:ion|atom)?$", re.I)
_COFACTOR = re.compile(
    r"\bha?eme?\b|porphyrin|chlorophyll|chlorin\b|pheophytin|caroten|lutein|xanthin|lycopene|spheroiden|rhodopin|"
    r"neurosporene|quinone|quinol\b|flavin|nicotinamide|retinal\b|retinol|cobalamin|pyridoxal|thiamin|molybd|cluster|"
    r"phycocyanobilin|phycoerythrobilin|phycourobilin|biliverdin|phytochromobilin|echinenone|spirilloxanthin|"
    r"lipoic|tocopher|plastoquin|\bcoenzyme q|bacteriochlor|dinuclear copper", re.I)
_NUCLEOTIDE = re.compile(r"adenosine|guanosine|uridine|cytidine|thymidine|inosine|adenyl|guanyl|uridyl|cytidyl|"
                         r"\bxanthosine|aminophosphonic acid", re.I)
_ALK = r"(?:hex|hept|oct|non|dec|undec|dodec|tridec|tetradec|hexadec|octadec|laur)"
_DETERGENT = re.compile(
    rf"{_ALK}yl\b[^|]*(?:gluco|malto|galacto|allo|altro)\w*side|maltose neopentyl|glucose neopentyl|neopentyl glycol|"
    r"digitonin|diosgenin|spirost|\bchapso?\b|cholamidopropyl|ammonio-1-propanesulfonate|fos-?choline|triton|"
    r"\btween\b|zwittergent|cymal|\bdeoxycholate\b|\bcholate\b|cholic acid|dimethylamine[- ]n-oxide|"
    r"(?:lauryl|dodecyl|decyl) sulfate|glycol mono\w*yl ether|polyoxyethylene|^\(hydroxyethyloxy\)|"
    rf"o-{_ALK}\w* \w*ethylene glycol|\bhega-\d|\bmega-\d|pentaoxatricosan|\bdetergent\b", re.I)
_STEROL = re.compile(r"cholest|sterol\b|ergost|stigmast|hemisuccinate", re.I)
_LIPID = re.compile(
    r"\b\d{1,2}:\d\b|phosphatidyl|phospho(?:choline|ethanolamine|serine|glycerol|inositol)|phosphatidic|"
    r"glycero-?3|sn-glycer|glycerol?\b.*(?:oyl|ate|acyl|anyl)\b|(?:oyl|acyl|anyl)\b.*glycerol|"
    r"oleoyl|palmitoyl|stearoyl|myristoyl|lauroyl|linoleoyl|arachidonoyl|"
    r"cardiolipin|ceramide|sphingo|ganglioside|diacylglycerol|diglycerid|\blipid\b|lipid [a-z]\b|"
    r"(?:^|\| )n?-?(?:pent|hex|hept|oct|non|dec|undec|dodec|tridec|tetradec|pentadec|hexadec|heptadec|octadec|eicos|"
    r"icos|docos)(?:ane|ene|an-1-ol|anol|anal)(?:$| \|)|"
    r"(?:octan|nonan|decan|undecan|dodecan|tridecan|tetradecan|pentadecan|hexadecan|heptadecan|octadecan|"
    r"octadec\w*en|hexadec\w*en|icos\w*|eicos\w*|docos\w*|cosan)(?:oyl|oate|oic acid)|"
    r"(?:oleic|stearic|palmitic|palmitoleic|myristic|lauric|linoleic|linolenic|arachidonic|petroselinic|elaidic|"
    r"vaccenic|capric|caprylic) acid|\b(?:ole|stear|palmit|myrist|laur)ate\b|galactosyl ?diacyl|"
    r"squalene|undecaprenyl|dolich|farnesyl|geranyl|phytanyl", re.I)
_ADDITIVE = re.compile(
    r"^glycerol$|^glycerin|ethanediol$|^ethane-1,2-diol|ethylene glycol$|hydroxyethyl\)ether|^sulfate ion|"
    r"^phosphate ion|\bacetate ion|^acetic acid|^acetyl group|"
    r"^formic acid|^formate|^nitrate|^citrate|^citric acid|tartrate$|tartaric acid$|^malonate|^malonic|^tris\b|"
    r"^2-amino-2-hydroxymethyl-propane-1,3-diol|hydroxymethyl.aminomethane|piperazine ethanesulfonic|"
    r"morpholino\)-(?:ethane|propane)sulfonic|^imidazole$|"
    r"^dimethyl sulfoxide|mercaptoethanol$|dithiothreitol|dithiane-diol|^azide|^thiocyanate|^ethanol$|^methanol$|"
    r"^isopropyl alcohol|^propan-2-ol|^[\w,()-]*(?:propane|butane|pentane|hexane)diol$|^(?:bi)?carbonate|"
    r"^ammonium|^urea$|ethylenediaminetetra|^oxalate|^succinic acid|^acetone$|carboxyethyl.phosphine|"
    r"^cacodylate|^arsenate|^borate|^boric|^pyrophosphate|^diphosphate$|^perchlorate|^sulfite|"
    r"unknown (?:ligand|atom)|^water$|^hydroxide|^peroxide|^oxygen|^(?:carbon|sulfur|nitrogen) (?:di|mon)oxide|"
    r"^nitric oxide|^ethanolamine$|^guanidine$|^glycine betaine|^trimethylamine oxide", re.I)
# The same kind of molecule known by a trade name or abbreviation among its synonyms.
_ADDITIVE_SYN = re.compile(r"dithiothreitol|\bEDTA\b|\bEGTA\b|\bTCEP\b|\bHEPES\b|\bDMSO\b|polyethylene glycol|"
                           r"\bspermidine\b|\bputrescine\b", re.I)
_ALKYL = re.compile(r"(?:C[=/\\]{0,3}){8,}")           # an unbranched chain of 8+ carbons in a SMILES string


def _elements(formula: str | None) -> dict[str, int]:
    return {m[1]: int(m[2] or 1) for m in re.finditer(r"([A-Z][a-z]?)(\d*)", (formula or "").split(" +")[0])
            if m[1]}


def classify(comp: dict) -> str:
    """Class of a chemical component from its own record ({name, formula,
    weight, type, smiles, synonyms}; anything but the name may be missing):
    ion, cofactor, nucleotide, detergent, lipid, sugar, additive, amino acid
    or drug (any other small molecule)."""
    name = (comp.get("name") or "").strip()
    names = " | ".join([name] + [s[0] for s in comp.get("synonyms") or [] if s[1] in ("PDB Reference Data", "Author")])
    el = _elements(comp.get("formula"))
    smiles = re.sub(r"Cl|Br", "X", comp.get("smiles") or "")
    weight = comp.get("weight") or 0
    ctype = (comp.get("type") or "").lower()
    organic = not el or "C" in el
    if (el and sum(el.values()) == 1) or (not el and _ELEMENT_ION.match(name)):
        return "ion"
    if _COFACTOR.search(names):
        return "cofactor"
    if _DETERGENT.search(names):
        return "detergent"
    if _NUCLEOTIDE.search(names):
        return "nucleotide"
    if _STEROL.search(names) or _LIPID.search(names):
        return "lipid"
    if "saccharide" in ctype:
        return "sugar"
    if set(el) & _METALS and (el.get("C", 0) >= 4 or (not organic and set(el) & {"Fe", "Cu", "Mo", "Ni", "Mn", "Co", "W"})):
        return "cofactor"  # a metal held in an organic ring or a metal cluster
    polyol = set(el) <= {"C", "H", "O"} and el.get("O", 0) >= 3 and el.get("C", 99) <= 2.5 * el["O"] and \
        smiles and not re.search(r"[=\d]", smiles)   # glycerol, PEG fragments
    if _ADDITIVE.search(name) or _ADDITIVE_SYN.search(names) or (el and not organic) or polyol:
        return "additive"
    if smiles and _ALKYL.search(smiles):
        return "lipid"
    if "peptide linking" in ctype:
        return "amino acid"
    if weight and weight < 90:
        return "additive"
    return "drug"


# ---------------------------------------------------------------- does a text mention a ligand?

_GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "Δ": "delta", "ε": "epsilon", "κ": "kappa",
          "μ": "mu", "ω": "omega", "λ": "lambda"}
# Words a synonym list may hold that name no particular molecule.
_GENERIC = set("ion ions acid water atom unknown ligand fragment group peptide lipid detergent sugar inhibitor "
               "agonist antagonist compound free bound form oxygen base salt ester ether amine alcohol".split())
# 'without ATP', 'in the absence of Ca2+', 'nucleotide-free', 'Ca2+-depleted': not a bound ligand.
_NEGATED = re.compile(r"\b(?:without|absence of|no|lacking|free of|devoid of|depleted of|minus)\s+(?:bound\s+|added\s+)?"
                      r"[\w+.-]+(?:\s*(?:,|/|and|or)\s*[\w+.-]+)*|[\w+.]+(?:-|\s)(?:free|depleted)\b", re.I)


def plain(s: str | None) -> str:
    """Markup-free text with Greek letters spelled out ('ATP&gamma;S' -> 'ATPgammaS')."""
    s = html.unescape(re.sub(r"<[^>]+>", "", s or ""))
    return "".join(_GREEK.get(ch, ch) for ch in s)


def _forms(s: str) -> tuple[str, str]:
    """Two comparable spellings of a text: words apart ('amp pnp') and with
    hyphens closed up ('amppnp'), each padded with a space at both ends."""
    s = re.sub(r"(?<![a-z0-9])\+", " ", plain(s).lower())   # '+' only as a charge ('ca2+'), not '+ GABA'
    apart = re.sub(r"[^a-z0-9+]+", " ", s).strip()
    closed = re.sub(r"[^a-z0-9+]+", " ", re.sub(r"(?<=[a-z0-9])[-‐–](?=[a-z0-9])", "", s)).strip()
    return f" {apart} ", f" {closed} "


def _usable(term: str) -> bool:
    """Long enough to mean one molecule: four letters or digits, or a three-character
    abbreviation or code ('ATP', 'IP3', '5-HT')."""
    t = plain(term).strip()
    core = re.sub(r"[^A-Za-z0-9]", "", t)
    if t.lower() in _GENERIC or not re.search(r"[A-Za-z]", core):
        return False
    return len(core) >= 4 or (len(core) == 3 and (sum(c.isupper() for c in core) >= 2 or any(c.isdigit() for c in core)))


def comp_terms(comp: dict, description: str = "") -> list[str]:
    """Every name a text may use for a chemical component: its PDB name and
    description, the PDB's and DrugBank's synonyms, the short name mpstruc
    uses (listing.COMMON) and plainer variants ('(S)-nicotine' -> 'nicotine',
    '... acid' -> '...ate')."""
    from . import listing
    cid = comp.get("id") or ""
    raw = [*listing.COMMON_ALT.get(cid, []), listing.COMMON.get(cid), description, comp.get("name"), comp.get("drug")]
    raw += [s[0] for s in comp.get("synonyms") or [] if s[1] in ("PDB Reference Data", "Author", "DrugBank")]
    out: list[str] = []
    for w in raw:
        for v in listing._variants(plain(w)) if w else []:
            if _usable(v) and v not in out:
                out.append(v)
    return out


def mentions(text: str, comp: dict, description: str = "", terms: list[str] | None = None) -> str | None:
    """The way `text` names the chemical component, or None if it does not.

    A name, synonym or abbreviation must stand as whole words (hyphens may
    differ: 'AMPPNP' = 'AMP-PNP'); the bare chemical ID counts only in capitals
    ('ATP', never 'act' for ACT). A single-atom ion counts when the text speaks
    of the ion ('Ca2+', 'zinc-bound'), not of a protein named after it
    ('calcium channel'). 'X-free' and 'without X' are not mentions."""
    from . import listing
    cid = comp.get("id") or ""
    text = _NEGATED.sub(" ", plain(text))
    if cid in listing.IONS:
        m = listing.ion_in_text(cid, text)
        return m[0] if m else None
    if len(cid) >= 3 and re.search(r"[A-Z]", cid) and cid not in ID_WORDS and \
            re.search(rf"(?<![A-Za-z0-9]){re.escape(cid)}(?![A-Za-z0-9])", text):
        return cid
    label = listing.COMMON.get(cid) or ""
    if re.fullmatch(r"[A-Z][a-z]{2}", label) and re.search(rf"(?<![A-Za-z0-9-]){label}(?![A-Za-z0-9-])", text):
        return label  # a free amino acid by its three-letter code ('with bound Arg', not 'Arg123' or 'Ala-Gln')
    forms = _forms(text)
    for t in terms if terms is not None else comp_terms(comp, description):
        if any(k in f for k in _forms(t) for f in forms):
            return t
    return None


# Chemical IDs that are also everyday abbreviations in mpstruc text (a positive allosteric
# modulator is a 'PAM', not palmitoleic acid), so the bare ID proves nothing.
ID_WORDS = set("PAM NAM TMD ECD LBD NTD CTD NBD FAB APO MSP WT1 ACT TOP DOM SUB MAP EM1 PLC PKA PKC ATD ICD VSD "
               "TRP ARG LYS HIS SER THR TYR PHE LEU ILE VAL MET PRO ALA GLY GLU ASP ASN GLN CYS SAM CAM DAG PIP".split())


# ---------------------------------------------------------------- what curators wrote for each structure

APO_TITLE = re.compile(r"(?<![A-Za-z-])(?:apo|unliganded|ligand[- ]free|substrate[- ]free|nucleotide[- ]free|"
                       r"empty|unbound|drug[- ]free|inhibitor[- ]free|agonist[- ]free)(?![A-Za-z])", re.I)
_BOUND_IN_NAME = re.compile(r"[\w+.()\[\]',&;<>/-]+[- ](?:bound|loaded|occupied|liganded)\b", re.I)


def curator_texts(db) -> dict[str, dict]:
    """PDB code -> what curators wrote about that structure's condition:
    {entry, family, kind, text, own} (`own` is the text without what a related
    structure inherits from its entry). For an entry: the condition part of its name
    ('with bound ATP'; the protein's own name is left out, so 'ATP synthase'
    names no ligand, but 'ATP-bound' inside it does) and the description lines
    that are not about a related structure. For a related structure: its own
    line, plus the entry's condition, which the line leaves unsaid when shared."""
    from . import listing
    out: dict[str, dict] = {}
    for e in db.entries():
        name = re.sub(r"\s+", " ", e.name or "").strip()
        base = listing.base_name(name) or ""
        cond = name[len(base):] if name.startswith(base) else ""
        cond += " ; " + " ; ".join(m[0] for m in _BOUND_IN_NAME.finditer(base))
        lines = [l.strip() for l in re.split(r"<br\s*/?>", e.description or "") if l.strip()]
        codes = {r.pdb.upper() for r in e.related}
        own = []
        rel: dict[str, str] = {}
        for l in lines:
            hit = [c for c in re.findall(r"(?<![A-Za-z0-9])[0-9][A-Za-z0-9]{3}(?![A-Za-z0-9])", l) if c.upper() in codes]
            if hit:
                for c in hit:
                    rel[c.upper()] = re.sub(r"[,;]?\s*[\d.]+\s*(?:Å|&Aring;|A)\b.*$", "", l)
            else:
                own.append(l)
        family = (e.masterProteinPdbCode or e.pdbCode).upper() if e.kind == "daughter" else e.pdbCode.upper()
        text = cond + " ; " + " ; ".join(own)
        out[e.pdbCode.upper()] = {"entry": e.pdbCode.upper(), "family": family, "kind": "entry", "own": text,
                                  "text": text, "year": e.bibliography.year}
        for c in codes:
            if c in out and out[c]["kind"] == "entry":
                continue
            out[c] = {"entry": e.pdbCode.upper(), "family": family, "kind": "related", "own": rel.get(c, ""),
                      "text": rel.get(c, "") + " ; " + cond, "year": e.bibliography.year}
    return out


def _comp_entity(comp: dict, description: str = "") -> dict:
    """A component as the non-polymer entity record listing.py reads."""
    return {"rcsb_nonpolymer_entity": {"pdbx_description": description},
            "nonpolymer_comp": {"chem_comp": {"id": comp.get("id"), "name": comp.get("name")}, "record": comp,
                                "rcsb_chem_comp_synonyms": [{"name": s[0], "provenance_source": s[1]}
                                                            for s in comp.get("synonyms") or []]}}


def observe(db, corpus: dict[str, dict | None] | None = None, comps: dict[str, dict | None] | None = None) -> list[dict]:
    """One record per (structure in mpstruc, non-polymer ligand in its PDB entry):
    {pdb, entry, family, kind, id, cls, named, how, title, cit, apo}.

    `named`: the curators' text for the structure mentions the ligand. `how`
    says on what evidence: 'direct' (its ID, name or a synonym is in the text),
    'fit' (the text names a compound by a code or abbreviation that fits this
    ligand's chemical name) or 'sole' (the text names one more compound and
    this is the structure's only other small-molecule ligand). `title` / `cit`:
    the PDB title / the paper's title mentions it; `apo`: the title says the
    structure is unliganded."""
    from . import listing
    texts = curator_texts(db)
    corpus = corpus if corpus is not None else corpus_entries(sorted(texts))
    ids = sorted({((n.get("nonpolymer_comp") or {}).get("chem_comp") or {}).get("id") for v in corpus.values() if v
                  for n in v.get("nonpolymer_entities") or []} - {None})
    comps = comps if comps is not None else chem_comps(ids)
    cls_of = {c: classify(comps[c] or {"id": c}) for c in ids}
    out = []
    for pdb, ctx in sorted(texts.items()):
        ent = corpus.get(pdb)
        if not ent:
            continue
        title = (ent.get("struct") or {}).get("title") or ""
        cit = (ent.get("rcsb_primary_citation") or {}).get("title") or ""
        rows, seen = [], set()
        for n in ent.get("nonpolymer_entities") or []:
            cid = ((n.get("nonpolymer_comp") or {}).get("chem_comp") or {}).get("id")
            if not cid or cid in HARD_EXCLUDE or cid in seen:
                continue
            seen.add(cid)
            comp = comps.get(cid) or {"id": cid, "name": ((n.get("nonpolymer_comp") or {}).get("chem_comp") or {}).get("name")}
            desc = (n.get("rcsb_nonpolymer_entity") or {}).get("pdbx_description") or ""
            terms = comp_terms(comp, desc)
            hit = mentions(ctx["text"], comp, desc, terms)
            rows.append({"pdb": pdb, "entry": ctx["entry"], "family": ctx["family"], "kind": ctx["kind"], "id": cid,
                         "cls": cls_of[cid], "named": bool(hit), "how": "direct" if hit else None, "by": hit,
                         "title": False,
                         "cit": bool(mentions(cit, comp, desc, terms)), "apo": bool(APO_TITLE.search(title)),
                         "_names": [comp.get("name") or "", desc] + [s[0] for s in comp.get("synonyms") or []],
                         "_w": comp.get("weight") or 0})
        in_title = listing.ligands_in_title({"polymer_entities": ent.get("polymer_entities"), "nonpolymer_entities": [
            _comp_entity(comps.get(r["id"]) or {"id": r["id"], "name": r["_names"][0]}, r["_names"][1]) for r in rows]}, title)
        for r in rows:
            r["title"] = r["id"] in in_title
        _assign_leftovers(rows, ctx["own"], ent, title)
        for r in rows:
            r.pop("_names"), r.pop("_w")
        out += rows
    return out


# Lipid words of curators' text and what the ligand's chemical names must then contain.
_LIPID_WORDS = [
    (r"(?:PO|DO|DM|DP|DS|DL|SO|SA)?PC|phosphatidylcholine", r"choline|trimethylammonio|trimethylazaniumyl"),
    (r"(?:PO|DO|DM|DP|DS|DL|SO|SA)?PE|phosphatidylethanolamine", r"ethanolamine|aminoethoxy|aminoethyl|azanylethoxy"),
    (r"(?:PO|DO|DM|DP|DS|DL|SO|SA)?PS|phosphatidylserine", r"serine"),
    (r"(?:PO|DO|DM|DP|DS|DL|SO|SA)?PG|phosphatidylglycerol", r"glycerol|dihydroxypropyl|bis\(oxidanyl\)propoxy"),
    (r"PIP\d?|PI|PI\(?[\d,]+\)?P\d?|PtdIns\S*|PI\d?P|phosphatidylinositol\S*|phosphoinositide", r"inositol|cyclohexyl"),
    (r"cholesterol|CHS|sterol", r"cholest|sterol"),
    (r"cardiolipin|CL|CDL", r"cardiolipin"),
    (r"LP[ACES]|S1P|DAG|diacylglycerol|(?:phospho|sphingo|glyco)?lipids?|sphingomyelin|ceramide|fatty acids?", r"."),
]
_ALKANE = re.compile(r"^n?-?[a-z]+ane$", re.I)


_BINDING = re.compile(r"(?:with bound|bound to|bound with|in complex with|complexed with|inhibited by|activated by|"
                      r"in (?:the )?presence of)\s+([^;,()]+)", re.I)


def _bound_phrases(text: str) -> list[str]:
    """Compounds and partners a curator's text names as bound: listing.title_names,
    plus the last word of a longer phrase when it reads as a code or name
    ('with bound a full bitopic agonist C5-Quino' -> 'C5-Quino') and
    'TC-I 2014-bound' -> 'TC-I 2014'."""
    from . import listing
    out = listing.title_names(text)
    extra = [m[1] for m in re.finditer(r"(?<![\w-])([A-Z][\w.+-]*(?: \d+)?)-bound\b", text)]
    for m in _BINDING.finditer(text):
        for part in re.split(r"\s+and\s+|\s*/\s*|\s+plus\s+|\s*\+\s*", m[1]):
            words = part.strip(" .").split()
            if len(words) >= 3 and re.search(r"[A-Z0-9]", words[-1][1:] or "") and not listing.PROTEINISH.search(words[-1]):
                extra.append(words[-1])
    for x in extra:
        if not any(x in o or o in x for o in out):
            out.append(x)
    return out


def _stems_fit(term: str, names: list[str]) -> bool:
    """A common name for one of the chemical names: two separate four-letter
    stems, covering at least half of it, are found in one name
    ('carbamylcholine' in 'carbamoylcholine')."""
    letters = re.sub(r"[^a-z]", "", term.lower())
    for n in (re.sub(r"[^a-z]", "", plain(x).lower()) for x in names):
        hits, i = 0, 0
        while i <= len(letters) - 4:
            if letters[i:i + 4] in n:
                hits, i = hits + 1, i + 4
            else:
                i += 1
        if hits >= 2 and 8 * hits >= len(letters):
            return True
    return False


def _assign_leftovers(rows: list[dict], text: str, ent: dict, title: str) -> None:
    """Compounds the curators' own text for the structure names that no
    ligand's names explain ('with bound LY354740' when the PDB only has a
    systematic name). In turn: a ligand whose chemical name shares stems with
    the word ('fit'); for a lipid word ('PIP2', 'POPC', 'lipid'), the structure's
    only lipid that is not a bare alkyl chain ('fit'); and when at least as
    many unexplained compounds are named as the structure has unexplained
    small-molecule ligands, those ligands; when fewer, as many ligands as there
    are compounds ('sole')."""
    from . import listing
    forms_named = [_forms(str(r["by"])) for r in rows if r["named"]]
    left = []
    for ph in _bound_phrases(plain(text)):
        k = _forms(ph)
        if any(k[0] in f[0] or f[0] in k[0] or k[1] in f[1] or f[1] in k[1] for f in forms_named):
            continue  # already explained by a ligand named directly
        if listing.names_protein(ph, ent) or listing._is_antibody(ph) or \
                listing.NOT_A_LIGAND.match(ph) or listing.MUTATION.fullmatch(ph):
            continue  # a protein partner or a word for the sample
        left.append(ph)
    codes = []
    for ph in left:
        free = [r for r in rows if not r["named"]]
        pick = next((r for r in free if r["cls"] in ("drug", "nucleotide", "amino acid", "lipid", "sugar", "cofactor")
                     and len(re.sub(r"[^a-z]", "", ph.lower())) >= 6 and _stems_fit(ph, r["_names"])), None)
        need = next((must for word, must in _LIPID_WORDS if re.fullmatch(word, ph, re.I)), None)
        if not pick and need:
            lipids = [r for r in free if r["cls"] == "lipid" and not _ALKANE.match(r["_names"][0])
                      and re.search(need, " ".join(r["_names"]), re.I)]
            pick = lipids[0] if len(lipids) == 1 else None
            if not pick:
                continue
        if pick:
            pick.update(named=True, how="fit", by=ph)
        elif len(ph.split()) <= 4 and not listing.SETTING.match(ph):
            codes.append(ph)
    # Small molecules first, then free amino acids and nucleotides; within a class the
    # ones the PDB title names, then the heaviest (a drug outweighs a stray additive).
    order = ("drug", "amino acid", "nucleotide")
    drugs = sorted((r for r in rows if not r["named"] and r["cls"] in order),
                   key=lambda r: (order.index(r["cls"]), not r["title"], -r["_w"]))[:len(codes)]
    # A free amino acid or nucleotide only for a compound the text clearly binds ('TNCA-bound').
    strong = [ph for ph in codes if re.search(
        rf"(?:bound|complex with|complexed with)\s+(?:\S+\s+){{0,3}}?{re.escape(ph)}|{re.escape(ph)}[- ]bound", text, re.I)]
    drugs = [r for r in drugs if r["cls"] == "drug" or strong]
    for r, ph in zip(drugs, codes + codes[-1:] * len(drugs)):
        r.update(named=True, how="sole", by=ph)


# ---------------------------------------------------------------- learned usage

def _tally(rows: list[dict]) -> dict:
    """n/named over all rows, and the same where the PDB title mentions the ligand."""
    t = [r for r in rows if r["title"]]
    return {"n": len(rows), "named": sum(r["named"] for r in rows),
            "n_title": len(t), "named_title": sum(r["named"] for r in t)}


RARE = 5  # a component seen in at most this many structures stands for 'one nobody has seen'


def summarise(obs: list[dict], comps: dict[str, dict | None] | None = None, rule: dict | None = None,
              meta: dict | None = None) -> dict:
    """Observations -> the ligand_usage.json document:
    comps    {id: {name, class, n, named, n_title, named_title}}
    classes  {class: {all: {...}, rare: {...}}}; 'rare' counts only components seen
             in at most RARE structures, the best guide to a component never seen
    overall  the same over every ligand
    rule     {threshold, prior}: see decide()."""
    comps = comps or {}
    by_comp, by_cls = defaultdict(list), defaultdict(list)
    for r in obs:
        by_comp[r["id"]].append(r)
        by_cls[r["cls"]].append(r)
    doc = {"meta": meta or {}, "rule": rule or {"threshold": 0.5, "prior": 2.0},
           "overall": _tally(obs), "classes": {}, "comps": {}}
    for c in CLASSES:
        rows = by_cls.get(c, [])
        doc["classes"][c] = {"label": CLASS_LABEL[c], "all": _tally(rows),
                             "rare": _tally([r for r in rows if len(by_comp[r["id"]]) <= RARE])}
    for cid in sorted(by_comp):
        rows = by_comp[cid]
        doc["comps"][cid] = {"name": (comps.get(cid) or {}).get("name") or "", "class": rows[0]["cls"], **_tally(rows)}
    return doc


def dump(doc: dict, path: Path = USAGE_PATH) -> None:
    """One component per line, so a rebuild shows up as a readable diff."""
    head = {k: v for k, v in doc.items() if k != "comps"}
    s = json.dumps(head, indent=1, ensure_ascii=False)[:-2] + ',\n "comps": {\n'
    s += ",\n".join(f'  {json.dumps(k)}: {json.dumps(v, ensure_ascii=False)}' for k, v in doc["comps"].items())
    path.write_text(s + "\n }\n}\n", encoding="utf-8")


_USAGE: dict[str, dict] = {}


def load_usage(path: Path | str | None = None) -> dict:
    """The learned usage (ligand_usage.json; MPSTRUC_LIGAND_USAGE names another file)."""
    import os
    p = str(path or os.environ.get("MPSTRUC_LIGAND_USAGE") or USAGE_PATH)
    if p not in _USAGE:
        _USAGE[p] = json.loads(Path(p).read_text(encoding="utf-8")) if Path(p).exists() else \
            {"rule": {"threshold": 0.5, "prior": 2.0}, "overall": {}, "classes": {}, "comps": {}}
    return _USAGE[p]


def _split(t: dict | None, said: bool) -> tuple[int, int]:
    """(structures, named) for the situation: the PDB title mentions the ligand, or does not."""
    t = t or {}
    if said:
        return t.get("n_title", 0), t.get("named_title", 0)
    return t.get("n", 0) - t.get("n_title", 0), t.get("named", 0) - t.get("named_title", 0)


def rate(usage: dict, cid: str, cls: str, said: bool) -> dict:
    """How often curators named this component, and its class, in this situation.
    `p` is the component's own rate pulled towards its class's rate by `prior`
    pseudo-structures: a component seen 200 times speaks for itself, one seen
    once or never goes by its class."""
    n, m = _split((usage.get("comps") or {}).get(cid), said)
    c = (usage.get("classes") or {}).get(cls) or {}
    cn, cm = _split(c.get("rare"), said)
    an, am = _split(c.get("all"), said)
    on, om = _split(usage.get("overall"), said)
    overall = om / on if on else (0.9 if said else 0.1)
    # The class's rate among rarely seen components, itself steadied by the whole class and all ligands.
    class_rate = (cm + 2 * ((am + 2 * overall) / (an + 2))) / (cn + 2)
    k = (usage.get("rule") or {}).get("prior", 2.0)
    return {"p": (m + k * class_rate) / (n + k), "n": n, "named": m, "class_rate": class_rate, "class_n": cn or an}


def class_of(comp: dict, usage: dict | None = None) -> str:
    """The class of a component: as learned for one mpstruc has seen, else from its record."""
    usage = usage if usage is not None else load_usage()
    return ((usage.get("comps") or {}).get(comp.get("id") or "") or {}).get("class") or classify(comp)


def decide(comp: dict, said: bool, usage: dict | None = None, cls: str | None = None) -> dict:
    """Would curators name this ligand? -> {id, name, class, decision, rate, why}.

    `comp` is the chemical component ({id, name, ...}), `said` whether the PDB
    title mentions it (listing.ligands_in_title). Named when the rate at which
    curators named it in this situation (see rate()) reaches the threshold
    chosen on tuning data; water and unknown atoms are never named."""
    usage = usage if usage is not None else load_usage()
    cid, name = comp.get("id") or "", comp.get("name") or ""
    known = (usage.get("comps") or {}).get(cid) or {}
    cls = cls or class_of(comp, usage)
    out = {"id": cid, "name": name or known.get("name") or cid, "class": cls}
    if cid in HARD_EXCLUDE:
        return {**out, "decision": "omitted", "rate": 0.0, "why": "water or an unidentified atom: never a named ligand"}
    r = rate(usage, cid, cls, said)
    named = r["p"] >= (usage.get("rule") or {}).get("threshold", 0.5)
    where = "the PDB title mentions it" if said else "the PDB title does not mention it"
    when = "when the PDB title mentions them" if said else "when the PDB title does not mention them"
    label = CLASS_LABEL.get(cls, cls)
    if r["n"] >= 5:
        why = (f"curators named {out['name']} in {r['named']} of the {r['n']} mpstruc structures that have it "
               f"and where {where}")
    else:
        seen = f"; mpstruc has this one in {r['n']} such structure{'s' * (r['n'] != 1)}, named in {r['named']}" \
            if r["n"] else "; mpstruc has not seen this one"
        why = (f"{label}: curators name such ligands in {round(100 * r['class_rate'])}% of mpstruc structures "
               f"{when}{seen}")
    return {**out, "decision": "named" if named else "omitted", "rate": round(r["p"], 3), "why": why}


# ---------------------------------------------------------------- build and measure

def _fold(family: str) -> str:
    """A fixed 60/20/20 split by entry family (a parent, its daughters and all their
    related structures stay together, so near-identical structures never straddle it)."""
    import hashlib
    h = int(hashlib.md5(family.encode()).hexdigest(), 16) % 10
    return "train" if h < 6 else "tune" if h < 8 else "test"


def score(rows: list[dict], usage: dict, comps: dict | None = None) -> dict:
    """Precision, recall and F1 of 'ligand named' as decide() would call it, against the curators."""
    comps = comps or {}
    tp = fp = fn = 0
    for r in rows:
        d = decide(comps.get(r["id"]) or {"id": r["id"]}, r["title"], usage, cls=r["cls"])
        pred = d["decision"] == "named" and not (r["apo"] and not r["title"])
        tp += pred and r["named"]
        fp += pred and not r["named"]
        fn += (not pred) and r["named"]
    p, rc = (tp / (tp + fp) if tp + fp else 0.0), (tp / (tp + fn) if tp + fn else 0.0)
    return {"precision": round(p, 3), "recall": round(rc, 3), "f1": round(2 * p * rc / (p + rc), 3) if p + rc else 0.0,
            "ligands": len(rows), "named": sum(r["named"] for r in rows)}


def tune(obs: list[dict], comps: dict | None = None) -> tuple[dict, dict]:
    """Choose the threshold and the prior strength on the tuning families
    (usage learned from the training families), then measure once on the
    held-out families. -> (rule, held-out score)."""
    folds = defaultdict(list)
    for r in obs:
        folds[_fold(r["family"])].append(r)
    tried = []
    for prior in (1.0, 2.0, 4.0):
        for threshold in (0.3, 0.4, 0.5, 0.6):
            rule = {"threshold": threshold, "prior": prior}
            tried.append((score(folds["tune"], summarise(folds["train"], comps, rule), comps)["f1"], rule))
    top = max(f1 for f1, _ in tried)
    # The scores form a plateau: among settings within 0.005 of the best, keep the plainest
    # (an even-odds threshold, a prior of two structures), so a rebuild does not flip between equals.
    best = min((x for x in tried if x[0] >= top - 0.005),
               key=lambda x: (abs(x[1]["threshold"] - 0.5), abs(x[1]["prior"] - 2.0)))
    rule = best[1]
    held = score(folds["test"], summarise(folds["train"] + folds["tune"], comps, rule), comps)
    return rule, {"tune_f1": best[0], **held, "families": len({r["family"] for r in folds["test"]})}


def build(out: Path = USAGE_PATH, exclude: set[str] = frozenset(), xml: str | None = None) -> dict:
    """Learn ligand usage from the mpstruc snapshot and write it to `out`.
    `exclude` lists entry PDB codes whose whole families are left out (to
    score drafts of those entries with a file that has never seen them)."""
    from .model import load
    from .report import source_xml
    db = load(xml or source_xml())
    obs = observe(db)
    if exclude:
        fams = {r["family"] for r in obs if r["pdb"] in exclude or r["entry"] in exclude}
        obs = [r for r in obs if r["family"] not in fams]
    ids = sorted({r["id"] for r in obs})
    comps = chem_comps(ids)
    rule, held = tune(obs, comps)
    doc = summarise(obs, comps, rule, {
        "source": db.meta.get("timeStamp") or "", "structures": len({r["pdb"] for r in obs}),
        "ligands": len(obs), "held_out": held,
        "note": "Built by `python -m mpstruc.ligands build` from the mpstruc snapshot: for each chemical component "
                "(and class), in how many mpstruc structures it is a ligand (n) and in how many of those the "
                "curators' text names it (named); n_title / named_title count the structures whose PDB title "
                "mentions it."})
    dump(doc, out)
    return doc


def _pct(a: int, b: int) -> str:
    return f"{a}/{b} ({100 * a / b:.0f}%)" if b else "0/0"


def report(doc: dict) -> str:
    """The class table of a usage document, as text."""
    rows = [f"{'class':12s} {'title mentions it':>22s} {'title does not':>22s}   (rarely seen components only)"]
    for c, v in doc["classes"].items():
        (tn, tm), (nn, nm) = _split(v["all"], True), _split(v["all"], False)
        (rtn, rtm), (rnn, rnm) = _split(v["rare"], True), _split(v["rare"], False)
        rows.append(f"{c:12s} {_pct(tm, tn):>22s} {_pct(nm, nn):>22s}   {_pct(rtm, rtn):>16s} {_pct(rnm, rnn):>16s}")
    h = doc["meta"].get("held_out") or {}
    rows.append(f"rule: named when the rate reaches {doc['rule']['threshold']} (prior {doc['rule']['prior']}); "
                f"held-out families: precision {h.get('precision')}, recall {h.get('recall')}, F1 {h.get('f1')} "
                f"on {h.get('ligands')} ligands of {h.get('families')} families")
    return "\n".join(rows)


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="learn ligand usage from mpstruc and write ligand_usage.json")
    b.add_argument("--out", default=str(USAGE_PATH))
    b.add_argument("--exclude", help="file of entry PDB codes whose families are left out")
    sub.add_parser("show", help="print the class table of ligand_usage.json")
    args = ap.parse_args()
    if args.cmd == "build":
        exclude = set(Path(args.exclude).read_text().split()) if args.exclude else set()
        doc = build(Path(args.out), {x.upper() for x in exclude})
        print(f"{doc['meta']['ligands']} ligands in {doc['meta']['structures']} mpstruc structures, "
              f"{len(doc['comps'])} chemical components -> {args.out}")
    else:
        doc = load_usage()
    print(report(doc))


if __name__ == "__main__":
    main()
