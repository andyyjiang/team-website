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
  structures, named in M', split by what the PDB title says of the ligand
  (situation()): it mentions it ('title': by a name, or by a family word such
  as 'PIP2' or a common name that fits its chemical names), it may mean it
  ('guess': an abbreviation, or a compound the title names that is credited
  to this ligand only because no other is left), or neither ('none'). The
  result is mpstruc/ligand_usage.json.
* DECIDE (`decide`). For a ligand of a new structure: the rate at which
  curators named this component (or, for a component mpstruc has rarely or
  never seen, its class) in the same situation, compared with a threshold
  chosen on tuning families. Every decision comes with a plain-English reason
  for the curators. What the title says is listing.ligands_in_title, the same
  test in learning and deciding. (The paper's title was tried as a second
  place to look; it did not help.)
* MEASURE (`tune`, `score`). Families of entries are split 60/20/20: usage is
  learned on the first part, the rule chosen on the second, and precision
  and recall of 'ligand named' reported on the third, which nothing saw.
  'Named' is what the text matcher found, not a curator's own record: about a
  quarter of the positive labels are inferred by elimination (how='sole'), and
  where curators copy a title word the label and the 'guess' situation agree by
  construction. ligand_usage.json therefore also holds the score on labels read
  in the text only, and the score for each situation: almost all of the gain
  over naming every drug-like ligand is made where the title does not mention
  the ligand, and that is also where the rule is weakest (about one in two).

A ligand's name inside a protein's name ('ATP synthase', 'IP3 receptor',
'fucoxanthin chlorophyll a/c-binding protein') is not a mention, in titles and
in curators' text alike.

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
    "polymer_entities { rcsb_polymer_entity { pdbx_description } entity_poly { rcsb_sample_sequence_length } } "
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
    (a component or a released entry never changes in the ways used here).
    An id is cached as None only when the PDB answered and does not know it. A
    failed request (network down, or an answer without data) caches nothing: those
    ids are left out of the result, which still holds everything already cached,
    and are asked for again on the next call."""
    from . import rcsb
    path = rcsb.CACHE_DIR / cache_name
    cache: dict = json.loads(path.read_text()) if path.exists() else {}
    todo = sorted({i.upper() for i in ids if i} - cache.keys())
    failed = None
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        try:
            data = rcsb._post_json(rcsb.GRAPHQL, {"query": query, "variables": {"ids": chunk}})
            got = (data.get("data") or {}).get(field)
            if got is None:  # HTTP 200 with {errors: [...], data: null}: no answer, not 'unknown'
                raise RuntimeError(f"the PDB returned no data ({str(data.get('errors'))[:200]})")
        except Exception as err:  # nothing is cached for what was not fetched: it is asked for again next time
            failed = err
            break
        found = {key(e): slim(e) for e in got if e}
        for cid in chunk:
            cache[cid] = found.get(cid)  # None: the PDB answered and does not know this id
        if (i // batch) % 10 == 9 or i + batch >= len(todo):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
    if failed:
        if cache:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(cache))
        print(f"ligands: could not fetch part of {cache_name} ({failed}); using what is cached")
    return {i: cache[i.upper()] for i in ids if i and i.upper() in cache}


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
    """Title, citation title, protein names and chain lengths, and non-polymer ligands
    per PDB entry; cached in cache/ligand_corpus_v2.json."""
    return _cached_fetch("ligand_corpus_v2.json", ids, CORPUS_QUERY, "entries", lambda e: e.get("rcsb_id"), batch)


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
    r"neurosporene|quinone|quinol\b|(?<!thio)flavin|nicotinamide|retinal\b|retinol|cobalamin|pyridoxal|thiamin|molybd|cluster|"
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


# In text reduced by _forms(): what follows a ligand's name when it is part of a protein's
# name ('ATP synthase', 'IP3 receptor', 'fucoxanthin chlorophyll a/c-binding protein').
_PROTEIN_AFTER = re.compile(
    r"(?:(?:ions?|\d\w*|[a-z]) ){0,2}(?:synthases?|synthetases?|receptors?|channels?|transporters?|pumps?|exchangers?|"
    r"symporters?|antiporters?|cotransporters?|uniporters?|carriers?|importers?|exporters?|permeases?|translocators?|"
    r"[a-z0-9]*binding|gated|dependent|sensitive|sensing|sensors?|(?!(?:ph|b|rele|c|incre|dise)ases? )[a-z]{3,}ases?) |"
    r"(?:[a-z0-9]+ ){0,4}[a-z0-9]*binding (?:proteins?|cassette) ")


def _found(k: str, f: str) -> bool:
    """Is the name `k` in the text `f` (both as _forms() gives them) other than inside a protein's name?"""
    i = f.find(k)
    while i >= 0:
        if not _PROTEIN_AFTER.match(f, i + len(k)):
            return True
        i = f.find(k, i + 1)
    return False


def mentions(text: str, comp: dict, description: str = "", terms: list[str] | None = None) -> str | None:
    """The way `text` names the chemical component, or None if it does not.

    A name, synonym or abbreviation must stand as whole words (hyphens may
    differ: 'AMPPNP' = 'AMP-PNP'); the bare chemical ID counts only in capitals
    ('ATP', never 'act' for ACT). A single-atom ion counts when the text speaks
    of the ion ('Ca2+', 'zinc-bound'), not of a protein named after it
    ('calcium channel'); likewise any other name inside a protein's name ('ATP
    synthase', 'chlorophyll a/c-binding protein') is not the ligand. 'X-free'
    and 'without X' are not mentions."""
    from . import listing
    cid = comp.get("id") or ""
    text = _NEGATED.sub(" ", plain(text))
    if cid in listing.IONS:
        m = listing.ion_in_text(cid, text)
        return m[0] if m else None
    if len(cid) >= 3 and re.search(r"[A-Z]", cid) and cid not in ID_WORDS and \
            re.search(rf"(?-i:{listing.name_pattern(cid)})(?!{listing.NAME_IN_PROTEIN})", text):
        return cid
    if cid in listing.CASED_ALIAS and re.search(listing.CASED_ALIAS[cid], text):
        return re.search(listing.CASED_ALIAS[cid], text)[0]
    label = listing.COMMON.get(cid) or ""
    if re.fullmatch(r"[A-Z][a-z]{2}", label) and re.search(rf"(?<![A-Za-z0-9-]){label}(?![A-Za-z0-9-])", text):
        return label  # a free amino acid by its three-letter code ('with bound Arg', not 'Arg123' or 'Ala-Gln')
    forms = _forms(text)
    for t in terms if terms is not None else comp_terms(comp, description):
        if any(_found(k, f) for k in _forms(t) for f in forms):
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
    {pdb, entry, family, kind, id, cls, named, how, by, said, title, cit, apo}.

    `named`: the curators' text for the structure mentions the ligand. `how`
    says on what evidence: 'direct' (its ID, name or a synonym is in the text),
    'fit' (the text names a compound by a code or abbreviation that fits this
    ligand's chemical name) or 'sole' (the text names one more compound and
    this is the structure's only other small-molecule ligand). `said`: how the
    PDB title mentions it (situation(): 'title', 'guess' or 'none'); `title`:
    said == 'title'; `cit`: the paper's title mentions it; `apo`: the PDB title
    says the structure is unliganded."""
    from . import listing
    texts = curator_texts(db)
    corpus = corpus if corpus is not None else corpus_entries(sorted(texts))
    ids = sorted({((n.get("nonpolymer_comp") or {}).get("chem_comp") or {}).get("id") for v in corpus.values() if v
                  for n in v.get("nonpolymer_entities") or []} - {None})
    comps = comps if comps is not None else chem_comps(ids)
    cls_of = {c: classify(comps.get(c) or {"id": c}) for c in ids}
    n_structures: dict[str, int] = defaultdict(int)   # mpstruc structures with each component
    for pdb in texts:
        for cid in {((n.get("nonpolymer_comp") or {}).get("chem_comp") or {}).get("id")
                    for n in (corpus.get(pdb) or {}).get("nonpolymer_entities") or []}:
            n_structures[cid] += 1
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
                         "title": False, "said": "none",
                         "cit": bool(mentions(cit, comp, desc, terms)), "apo": bool(APO_TITLE.search(title)),
                         "_names": [comp.get("name") or "", desc] + [s[0] for s in comp.get("synonyms") or []],
                         "_w": comp.get("weight") or 0})
        in_title = listing.ligands_in_title({"polymer_entities": ent.get("polymer_entities"), "nonpolymer_entities": [
            _comp_entity(comps.get(r["id"]) or {"id": r["id"], "name": r["_names"][0]}, r["_names"][1]) for r in rows]},
            title, seen=n_structures, classes=cls_of)
        for r in rows:
            r["said"] = situation(in_title.get(r["id"]))   # 'title', 'guess' or 'none'
            r["title"] = r["said"] == "title"
        _assign_leftovers(rows, ctx["own"], ent, title, ctx["text"])
        for r in rows:
            r.pop("_names"), r.pop("_w")
        out += rows
    return out


# Words for a family of molecules ('PIP2', 'POPC', 'LPA', 'C22-CoA', 'ubiquinone') that a
# PDB title or a curator's text uses for a ligand the PDB knows by a systematic name:
# (the word, what the ligand's chemical names must then contain, the classes it may have).
_ACYL = r"(?:PO|DO|DM|DP|DS|DL|SO|SA|L|[Ll]yso-?|NBD-|\d\d:\d[- ])?"
_LIPIDS, _ANY = ("lipid",), ("lipid", "nucleotide", "drug", "cofactor")
FAMILY_WORDS = [
    (rf"(?-i:{_ACYL}PC)|PtdCho|(?:lyso-?)?phosphatidylcholines?", r"choline|trimethylammonio|trimethylazaniumyl", _LIPIDS),
    (rf"(?-i:{_ACYL}PE)|PtdEtn|(?:lyso-?)?phosphatidylethanolamines?", r"ethanolamine|aminoethoxy|aminoethyl|azanylethoxy", _LIPIDS),
    (rf"(?-i:{_ACYL}PS)|PtdSer|(?:lyso-?)?phosphatidylserines?", r"serine", _LIPIDS),
    (rf"(?-i:{_ACYL}PG)|PtdGro|(?:lyso-?)?phosphatidylglycerols?", r"glycerol|dihydroxypropyl|bis\(oxidanyl\)propoxy", _LIPIDS),
    (rf"(?-i:{_ACYL}PA)|(?:lyso-?)?phosphatidic acids?|phosphatidates?",
     r"^(?!.*(?:choline|ammonio|azaniumyl|ethanolamine|aminoethoxy|aminoethyl|serine|inositol)).*phosph", _LIPIDS),
    (r"PIP\d?|PI|L?PI\(?[\d,]+\)?P\d?|PtdIns\S*|PI\d?P|LPI|(?:lyso-?)?phosphatidylinositol\S*|phosphoinositides?",
     r"inositol|cyclohexyl", _LIPIDS),
    (r"cholesterol|CHS|sterol", r"cholest|sterol", _LIPIDS),
    (r"cardiolipin|CL|CDL", r"cardiolipin", _LIPIDS),
    (r"(?-i:[\w:,()-]*CoA)|[\w:,()-]*coenzyme a", r"coenzyme a|coa\b", _ANY),
    (r"[\w-]*(?:oyl|acyl)-(?:AMP|adenylate)", r"aminopurin|adenos|adenyl", _ANY),
    (r"\w*quino(?:ne|l)s?(?:[- ]?\d+)?|U?Q\d{1,2}|UQ|CoQ\d*|vitamin K\d?|M?K-?\d", r"quinon|quinol|naphthalenedione|naphthoquin",
     ("cofactor", "lipid", "drug")),
    (r"LP[ACES]|S1P|sphingosine[- ]1[- ]phosphate|DAG|diacylglycerol|(?:phospho|sphingo|glyco)?lipids?|sphingomyelin|"
     r"ceramide|fatty acids?", r".", _LIPIDS),
]
_ALKANE = re.compile(r"^n?-?[a-z]+ane$", re.I)


def _family_core(phrase: str) -> str:
    """'d16:1 S1P' -> 'S1P', 'L-alpha-lysophosphatidylinositol' -> 'lysophosphatidylinositol'."""
    return re.sub(r"^(?:[A-Za-z]?\d\d:\d(?:/\d\d:\d)?[- ]|(?:L-)?alpha-)+", "", plain(phrase).strip())


def is_family_word(phrase: str) -> bool:
    return any(re.fullmatch(word, _family_core(phrase), re.I) for word, _, _ in FAMILY_WORDS)


def family_picks(phrases: list[str], candidates: list[tuple]) -> dict:
    """{ligand key: phrase} for the family words among `phrases`: each family's
    words go to the structure's candidates (key, class, [chemical names]) of
    that family that are not bare alkyl chains, when there are no more such
    ligands than words for them ('C22-CoA and MCoA' with two acyl-CoA ligands;
    'PIP2' with one inositol lipid, never with two)."""
    out: dict = {}
    used: set[str] = set()
    for word, need, classes in FAMILY_WORDS:
        said = [ph for ph in dict.fromkeys(phrases) if ph not in used and re.fullmatch(word, _family_core(ph), re.I)]
        hits = [c for c in candidates if c[0] not in out and c[1] in classes
                and not _ALKANE.match(c[2][0] if c[2] else "") and re.search(need, " | ".join(x or "" for x in c[2]), re.I)]
        if hits and len(hits) <= len(said):
            for c, ph in zip(hits, said):
                out[c[0]] = ph
            used.update(said)
    return out


# Classes a common name from a text may fit by its stems (never an ion, a buffer or a detergent).
FIT_CLASSES = ("drug", "nucleotide", "amino acid", "lipid", "sugar", "cofactor")


_BINDING = re.compile(r"(?:with bound|bound to|bound with|in complex with|complexed with|inhibited by|activated by|"
                      r"in (?:the )?presence of)\s+([^;,()]+)", re.I)


# A function word ('With') or a phrase about protein chains ('5 Lhcr antennae', 'auxiliary subunit TARP gamma8').
_NOT_A_COMPOUND = re.compile(r"^(?:with|and|the|to|of|in|by|from|for|or|bound)$|"
                             r"\b(?:antennae?|subunits?|proteins?|chains?|domains?|protomers?|monomers?)\b", re.I)


def _bound_phrases(text: str) -> list[str]:
    """Compounds and partners a curator's text names as bound: listing.title_names,
    plus the last word of a longer phrase when it reads as a code or name
    ('with bound a full bitopic agonist C5-Quino' -> 'C5-Quino') and
    'TC-I 2014-bound' -> 'TC-I 2014'."""
    from . import listing
    out = listing.compound_phrases(text)
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


def stems_fit(term: str, names: list[str]) -> bool:
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


def _assign_leftovers(rows: list[dict], own: str, ent: dict, title: str, text: str | None = None) -> None:
    """Compounds the curators' text for the structure names that no ligand's
    names explain ('with bound LY354740' when the PDB only has a systematic
    name). In turn, how='fit': a ligand whose chemical name shares stems with
    the word; for a family word ('PIP2', 'POPC', 'LPA', 'acyl-CoA', 'lipid'), the
    structure's only ligand of that family that is not a bare alkyl chain; a
    ligand whose common name the word abbreviates ('NQNO'). These read the whole
    `text`. Then how='sole', from the structure's `own` text only (not what a
    related structure inherits from its entry): when at least as many
    unexplained compounds are named as the structure has unexplained
    small-molecule ligands, those ligands; when fewer, as many ligands as there
    are compounds."""
    from . import listing
    own_phrases = {_forms(ph)[0] for ph in _bound_phrases(plain(own))}
    left = []
    for ph in _bound_phrases(plain(text if text is not None else own)):
        k = _forms(ph)
        if any(k[0] in f[0] or f[0] in k[0] or k[1] in f[1] or f[1] in k[1]
               for f in (_forms(str(r["by"])) for r in rows if r["named"])):
            continue  # already explained by a ligand named directly
        if listing.names_protein(ph, ent) or listing._is_antibody(ph) or _NOT_A_COMPOUND.search(ph) or \
                listing.NOT_A_LIGAND.match(ph) or listing.MUTATION.fullmatch(ph) or listing.ION_NAME.fullmatch(ph.strip()):
            continue  # a protein partner, an ion or a word for the sample
        pick = next((r for r in rows if not r["named"] and r["cls"] in FIT_CLASSES
                     and len(re.sub(r"[^a-z]", "", ph.lower())) >= 6 and stems_fit(ph, r["_names"])), None)
        if pick:
            pick.update(named=True, how="fit", by=ph)
        else:
            left.append(ph)
    free = [r for r in rows if not r["named"]]
    for i, ph in family_picks(left, [(i, r["cls"], r["_names"]) for i, r in enumerate(free)]).items():
        free[i].update(named=True, how="fit", by=ph)
    codes = []
    for ph in left:
        if is_family_word(ph):
            continue   # a family word the ligands do not settle is nobody's
        pick = None
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9-]{1,6}", ph) and sum(c.isupper() for c in ph) >= 2:
            fits = [(listing._abbreviation_fit(ph, r["_names"]), r) for r in rows if not r["named"] and r["cls"] != "ion"]
            pick = max((x for x in fits if x[0]), key=lambda x: x[0], default=(0, None))[1]
        if pick:
            pick.update(named=True, how="fit", by=ph)
        elif len(ph.split()) <= 4 and not listing.SETTING.match(ph) and _forms(ph)[0] in own_phrases:
            codes.append(ph)
    # Small molecules first, then free amino acids and nucleotides; within a class the
    # ones the PDB title names, then the heaviest (a drug outweighs a stray additive).
    order = ("drug", "amino acid", "nucleotide")
    drugs = sorted((r for r in rows if not r["named"] and r["cls"] in order),
                   key=lambda r: (order.index(r["cls"]), r["said"] == "none", -r["_w"]))[:len(codes)]
    # A free amino acid or nucleotide only for a compound the text clearly binds ('TNCA-bound').
    strong = [ph for ph in codes if re.search(
        rf"(?:bound|complex with|complexed with)\s+(?:\S+\s+){{0,3}}?{re.escape(ph)}|{re.escape(ph)}[- ]bound", own, re.I)]
    drugs = [r for r in drugs if r["cls"] == "drug" or strong]
    for r, ph in zip(drugs, codes + codes[-1:] * len(drugs)):
        r.update(named=True, how="sole", by=ph)


# ---------------------------------------------------------------- learned usage

SITUATIONS = ("title", "guess", "none")
SITUATION_TEXT = {
    "title": ("the PDB title mentions it", "when the PDB title mentions them"),
    "guess": ("the PDB title names a compound that no ligand's name explains, which may be this one",
              "when the PDB title names a compound that may be them"),
    "none": ("the PDB title does not mention it", "when the PDB title does not mention them")}


def situation(how) -> str:
    """How firmly the PDB title mentions a ligand, from listing.ligands_in_title's
    answer: 'title' (by a name, its ID, or a family word or common name that fits
    its chemical names), 'guess' (an abbreviation that could be read from its
    name, or a compound the title names that is credited to it only because no
    other ligand is left), 'none'. The two are counted apart: curators name the
    ligand in most 'title' cases, far fewer 'guess' cases."""
    if how in (True, "name", "fit", "title"):
        return "title"
    return "guess" if how in ("abbreviation", "compound", "guess") else "none"


def _tally(rows: list[dict]) -> dict:
    """n/named over all rows, the same where the PDB title mentions the ligand (_title)
    and where the title only may mean it (_guess)."""
    t = [r for r in rows if r["said"] == "title"]
    g = [r for r in rows if r["said"] == "guess"]
    return {"n": len(rows), "named": sum(r["named"] for r in rows),
            "n_title": len(t), "named_title": sum(r["named"] for r in t),
            "n_guess": len(g), "named_guess": sum(r["named"] for r in g)}


RARE = 5  # a component seen in at most this many structures stands for 'one nobody has seen'


def summarise(obs: list[dict], comps: dict[str, dict | None] | None = None, rule: dict | None = None,
              meta: dict | None = None) -> dict:
    """Observations -> the ligand_usage.json document:
    comps    {id: {name, class, n, named, n_title, named_title, n_guess, named_guess}}
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


def _split(t: dict | None, said) -> tuple[int, int]:
    """(structures, named) for one situation: the PDB title mentions the ligand, may mean it, or neither."""
    t, sit = t or {}, situation(said)
    if sit == "title":
        return t.get("n_title", 0), t.get("named_title", 0)
    if sit == "guess":
        return t.get("n_guess", 0), t.get("named_guess", 0)
    return (t.get("n", 0) - t.get("n_title", 0) - t.get("n_guess", 0),
            t.get("named", 0) - t.get("named_title", 0) - t.get("named_guess", 0))


def rate(usage: dict, cid: str, cls: str, said) -> dict:
    """How often curators named this component, and its class, in this situation.
    `p` is the component's own rate pulled towards its class's rate by `prior`
    pseudo-structures: a component seen 200 times speaks for itself, one seen
    once or never goes by its class."""
    n, m = _split((usage.get("comps") or {}).get(cid), said)
    c = (usage.get("classes") or {}).get(cls) or {}
    cn, cm = _split(c.get("rare"), said)
    an, am = _split(c.get("all"), said)
    on, om = _split(usage.get("overall"), said)
    overall = om / on if on else {"title": 0.9, "guess": 0.5, "none": 0.1}[situation(said)]
    # The class's rate among rarely seen components, itself steadied by the whole class and all ligands.
    class_rate = (cm + 2 * ((am + 2 * overall) / (an + 2))) / (cn + 2)
    k = (usage.get("rule") or {}).get("prior", 2.0)
    return {"p": (m + k * class_rate) / (n + k), "n": n, "named": m, "class_rate": class_rate, "class_n": cn or an}


def class_of(comp: dict, usage: dict | None = None) -> str:
    """The class of a component: as learned for one mpstruc has seen, else from its record."""
    usage = usage if usage is not None else load_usage()
    return ((usage.get("comps") or {}).get(comp.get("id") or "") or {}).get("class") or classify(comp)


def decide(comp: dict, said, usage: dict | None = None, cls: str | None = None) -> dict:
    """Would curators name this ligand? -> {id, name, class, decision, rate, why}.

    `comp` is the chemical component ({id, name, ...}), `said` how the PDB
    title mentions it (listing.ligands_in_title's answer, see situation();
    True / False also do). Named when the rate at which
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
    where, when = SITUATION_TEXT[situation(said)]
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


def predict(r: dict, usage: dict, comps: dict | None = None) -> bool:
    """Would the rule name the ligand of observation `r`? (As listing.ligands does: an
    apo / unliganded title names nothing it does not itself mention.)"""
    d = decide((comps or {}).get(r["id"]) or {"id": r["id"]}, r["said"], usage, cls=r["cls"])
    return d["decision"] == "named" and not (r["apo"] and r["said"] == "none")


def score(rows: list[dict], usage: dict, comps: dict | None = None) -> dict:
    """Precision, recall and F1 of 'ligand named' as decide() would call it, against the curators."""
    comps = comps or {}
    tp = fp = fn = 0
    for r in rows:
        pred = predict(r, usage, comps)
        tp += pred and r["named"]
        fp += pred and not r["named"]
        fn += (not pred) and r["named"]
    p, rc = (tp / (tp + fp) if tp + fp else 0.0), (tp / (tp + fn) if tp + fn else 0.0)
    return {"precision": round(p, 3), "recall": round(rc, 3), "f1": round(2 * p * rc / (p + rc), 3) if p + rc else 0.0,
            "ligands": len(rows), "named": sum(r["named"] for r in rows)}


# Where the choice is open, a ligand is named when curators named comparable ones in two of
# five structures: a ligand named in error is struck out at a glance (ligandNotes says why it
# is there), a missing one has to be looked up. At 0.3 a single mpstruc structure that names
# a detergent is enough to name it everywhere; at 0.5 no drug the title does not mention is named.
LEAN = 0.4


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
            sc = score(folds["tune"], summarise(folds["train"], comps, rule), comps)
            tried.append((sc["f1"], sc["recall"], rule))
    top = max(x[0] for x in tried)
    # The scores form a plateau (F1 within 0.005 from a threshold of 0.3 to 0.6): what moves
    # along it is whether a drug-like ligand the title does not mention is named, which curators
    # do about four times in ten, and no split of the data settles that. So among equals the
    # threshold nearest LEAN is kept, then the plainest prior, and a rebuild does not flip
    # between equals.
    best = min((x for x in tried if x[0] >= top - 0.005),
               key=lambda x: (abs(x[2]["threshold"] - LEAN), abs(x[2]["prior"] - 2.0)))
    best = (best[0], best[2])
    rule = best[1]
    usage = summarise(folds["train"] + folds["tune"], comps, rule)
    test = folds["test"]
    held = score(test, usage, comps)
    brief = lambda rows: {k: v for k, v in score(rows, usage, comps).items() if k != "f1" or rows}
    return rule, {
        "tune_f1": best[0], **held, "families": len({r["family"] for r in test}),
        # The same held-out ligands by what the PDB title says of them: the gain over naming
        # nothing is made where the title does not mention the ligand, and is smallest there.
        "by_situation": {s: brief([r for r in test if r["said"] == s]) for s in SITUATIONS},
        # 'sole' labels are inferred by elimination, not read in the curators' text: the
        # score without them (learned and measured on text matches only) is the firmer one.
        "text_matches_only": score([r for r in test if r["how"] != "sole"],
                                   summarise([r for r in folds["train"] + folds["tune"] if r["how"] != "sole"],
                                             comps, rule), comps)}


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
        "labels": {"named": sum(r["named"] for r in obs),
                   **{f"by_{h}": sum(r["how"] == h for r in obs) for h in ("direct", "fit", "sole")}},
        "note": "Built by `python -m mpstruc.ligands build` from the mpstruc snapshot: for each chemical component "
                "(and class), in how many mpstruc structures it is a ligand (n) and in how many of those the "
                "curators' text names it (named); n_title / named_title count the structures whose PDB title "
                "mentions it by a name, n_guess / named_guess those where the title only names an abbreviation or "
                "an unexplained compound that may be it. 'named' is what a text matcher found, not a curator's "
                "own record: labels.by_direct were read in the curators' text, by_fit matched by a code or family "
                "word, by_sole inferred because the text names one more compound than the ligands explain. "
                "held_out.text_matches_only is the score without the inferred labels."})
    dump(doc, out)
    return doc


def _pct(a: int, b: int) -> str:
    return f"{a}/{b} ({100 * a / b:.0f}%)" if b else "0/0"


def report(doc: dict) -> str:
    """The class table of a usage document, as text."""
    rows = [f"{'class':11s} {'title mentions it':>17s} {'title may mean it':>17s} {'title does not':>17s}   "
            f"(rarely seen components only)"]
    for c, v in doc["classes"].items():
        cells = [_pct(*reversed(_split(v[k], s))) for k in ("all", "rare") for s in SITUATIONS]
        rows.append(f"{c:11s} {cells[0]:>17s} {cells[1]:>17s} {cells[2]:>17s}   {cells[3]:>15s} {cells[4]:>13s} {cells[5]:>14s}")
    h = doc["meta"].get("held_out") or {}
    rows.append(f"rule: named when the rate reaches {doc['rule']['threshold']} (prior {doc['rule']['prior']}); "
                f"held-out families: precision {h.get('precision')}, recall {h.get('recall')}, F1 {h.get('f1')} "
                f"on {h.get('ligands')} ligands of {h.get('families')} families")
    for s, x in (h.get("by_situation") or {}).items():
        rows.append(f"  {SITUATION_TEXT[s][0]}: precision {x.get('precision')}, recall {x.get('recall')}, "
                    f"{x.get('named')} named of {x.get('ligands')}")
    t = h.get("text_matches_only") or {}
    lab = doc["meta"].get("labels") or {}
    if t:
        rows.append(f"  labels read in the curators' text only ({lab.get('by_sole')} of {lab.get('named')} 'named' labels are "
                    f"inferred by elimination and left out): precision {t.get('precision')}, recall {t.get('recall')}, "
                    f"F1 {t.get('f1')}")
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
