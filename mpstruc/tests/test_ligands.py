"""Which ligands curators name: the matcher, the classes, the learned rule
(mpstruc/ligands.py), its use in drafted listings and `automation relist`."""
import json
from types import SimpleNamespace

import pytest

from mpstruc import automation, ligands
from mpstruc.listing import build, facets, ligands as listing_ligands, ligands_in_title

from test_listing import entity, entry, ligand

ATP = {"id": "ATP", "name": "ADENOSINE-5'-TRIPHOSPHATE", "synonyms": [["ATP", "DrugBank"]]}
ANP = {"id": "ANP", "name": "PHOSPHOAMINOPHOSPHONIC ACID-ADENYLATE ESTER"}
CA = {"id": "CA", "name": "CALCIUM ION", "formula": "Ca"}
NA = {"id": "NA", "name": "SODIUM ION", "formula": "Na"}
KR0 = {"id": "KR0", "name": "3-decanoyloxypropyl decanoate", "formula": "C23 H44 O4", "weight": 384.6,
       "type": "non-polymer", "smiles": "CCCCCCCCCC(=O)OCCCOC(=O)CCCCCCCCC"}
NICOTINE = {"id": "NCT", "name": "(S)-3-(1-METHYL-2-PYRROLIDINYL)PYRIDINE", "weight": 162.2,
            "synonyms": [["(S)-(-)-NICOTINE", "PDB Reference Data"], ["Nicotine", "DrugBank"]]}


def tally(n, named, n_title, named_title, n_guess=0, named_guess=0):
    return {"n": n, "named": named, "n_title": n_title, "named_title": named_title,
            "n_guess": n_guess, "named_guess": named_guess}


# A small usage table with the shape of ligand_usage.json.
USAGE = {
    "rule": {"threshold": 0.5, "prior": 2.0},
    "overall": tally(1000, 150, 110, 90),
    "classes": {
        "ion": {"all": tally(500, 50, 60, 45), "rare": tally(40, 8, 8, 8)},
        "lipid": {"all": tally(700, 20, 12, 11), "rare": tally(50, 7, 4, 4)},
        "drug": {"all": tally(300, 220, 120, 114, 60, 50), "rare": tally(200, 150, 70, 68, 40, 34)},
        "nucleotide": {"all": tally(120, 50, 55, 42, 20, 6), "rare": tally(10, 4, 3, 3, 2, 1)},
        "cofactor": {"all": tally(400, 12, 12, 9, 5, 1), "rare": tally(30, 4, 4, 3, 2, 0)},
        "sugar": {"all": tally(120, 4, 6, 4, 3, 1), "rare": tally(10, 2, 3, 2, 1, 0)},
        "detergent": {"all": tally(140, 1, 1, 1), "rare": tally(6, 0, 0, 0)},
        "additive": {"all": tally(200, 9, 14, 9), "rare": tally(16, 2, 3, 2)},
    },
    "comps": {
        "CA": {"name": "CALCIUM ION", "class": "ion", **tally(800, 200, 190, 168)},
        "NA": {"name": "SODIUM ION", "class": "ion", **tally(800, 55, 80, 50)},
        "CLR": {"name": "CHOLESTEROL", "class": "lipid", **tally(760, 23, 37, 20)},
        "ATP": {"name": "ADENOSINE-5'-TRIPHOSPHATE", "class": "nucleotide", **tally(415, 188, 301, 172)},
    },
}


@pytest.fixture
def usage(monkeypatch):
    monkeypatch.setattr(ligands, "load_usage", lambda path=None: USAGE)
    return USAGE


# ---------------------------------------------------------------- the matcher

def test_names_synonyms_and_spellings_match_as_whole_words():
    assert ligands.mentions("with bound AMP-PNP, inward-open", ANP)       # mpstruc's short name
    assert ligands.mentions("AMPPNP-bound state", ANP)                    # hyphens may differ
    assert ligands.mentions("in complex with ATP&gamma;S and nicotine", NICOTINE)
    assert ligands.mentions("+ GABA + Zn<sup>2+</sup>", {"id": "ABU", "name": "GAMMA-AMINO-BUTANOIC ACID",
                                                          "synonyms": [["GABA", "DrugBank"]]})
    assert not ligands.mentions("with bound nicotinamide", NICOTINE)      # part of another word
    assert not ligands.mentions("apo form", NICOTINE)


def test_a_bare_chemical_id_counts_only_in_capitals():
    act = {"id": "ACT", "name": "ACETATE ION"}
    assert ligands.mentions("E1-ATP state", ATP) == "ATP"
    assert not ligands.mentions("these act together", act)
    pam = {"id": "PAM", "name": "PALMITOLEIC ACID"}
    assert not ligands.mentions("with bound PAM VU0364770", pam)          # a positive allosteric modulator
    assert ligands.mentions("with bound palmitoleate", pam)


def test_ions_count_when_bound_not_when_they_name_the_protein():
    assert ligands.mentions("with bound Ca<sup>2+</sup>", CA)
    assert ligands.mentions("zinc-bound, 2.9 &Aring;", {"id": "ZN", "name": "ZINC ION"})
    assert not ligands.mentions("Ca<sup>2+</sup>-activated K<sup>+</sup> channel", CA)
    assert not ligands.mentions("similar to the NhaA Na<sup>+</sup>/H<sup>+</sup> antiporter", NA)
    assert not ligands.mentions("potassium-chloride co-transporter", {"id": "K", "name": "POTASSIUM ION"})
    assert not ligands.mentions("two pore domain potassium ion channel", {"id": "K", "name": "POTASSIUM ION"})


def test_free_and_without_are_not_mentions():
    assert not ligands.mentions("nucleotide-free, ATP-free state", ATP)
    assert not ligands.mentions("in the absence of Ca<sup>2+</sup>", CA)
    assert not ligands.mentions("without ATP", ATP)
    assert ligands.mentions("with ATP, without ADP", ATP)


def test_free_amino_acids_by_three_letter_code():
    arg = {"id": "ARG", "name": "ARGININE"}
    assert ligands.mentions("in complex with 4F2hc, with bound Arg", arg) == "Arg"
    assert not ligands.mentions("Arg123 mutant", arg)
    assert not ligands.mentions("with Ala-Arg peptide", arg)


def db_of(*entries):
    return SimpleNamespace(entries=lambda: list(entries), meta={})


def mp_entry(pdb, name, description="", related=(), kind="parent", master=None):
    return SimpleNamespace(pdbCode=pdb, name=name, description=description, kind=kind, masterProteinPdbCode=master,
                           related=[SimpleNamespace(pdb=r) for r in related], bibliography=SimpleNamespace(year="2025"))


def corpus_entry(title, comps, proteins=("Receptor X",)):
    return {"struct": {"title": title}, "rcsb_primary_citation": {"title": ""},
            "polymer_entities": [{"rcsb_polymer_entity": {"pdbx_description": p}} for p in proteins],
            "nonpolymer_entities": [{"rcsb_nonpolymer_entity": {"pdbx_description": c["name"]},
                                     "nonpolymer_comp": {"chem_comp": {"id": c["id"], "name": c["name"]}}} for c in comps]}


def test_curator_text_is_the_condition_not_the_protein_name():
    db = db_of(mp_entry("9AAA", "ATP synthase, rotational state 1", "cryo-EM structure<br /> state 2, 3.1 &Aring;: 9AAB",
                        ["9AAB"]),
               mp_entry("9AAC", "P2X<sub>3</sub> receptor with bound ATP", "cryo-EM structure<br /> apo form, 3.0 &Aring;: 9AAD",
                        ["9AAD"], kind="daughter", master="9AAA"))
    t = ligands.curator_texts(db)
    assert not ligands.mentions(t["9AAA"]["text"], ATP)          # 'ATP synthase' is the protein
    assert ligands.mentions(t["9AAC"]["text"], ATP)
    assert ligands.mentions(t["9AAD"]["text"], ATP)              # a related line inherits its entry's condition ...
    assert not ligands.mentions(t["9AAD"]["own"], ATP)           # ... but its own words are just 'apo form'
    assert {t[p]["family"] for p in t} == {"9AAA"}              # parent, daughter and related structures: one family


def test_observe_a_named_drug_by_code_and_an_unnamed_lipid(usage):
    drug = {"id": "40F", "name": "(1S,2S,5R,6S)-2-aminobicyclo[3.1.0]hexane-2,6-dicarboxylic acid", "weight": 185.2}
    clr = {"id": "CLR", "name": "CHOLESTEROL"}
    db = db_of(mp_entry("9AAA", "mGlu<sub>2</sub> receptor with bound LY354740",
                        "cryo-EM structure<br /> apo form, 3.10 &Aring;: 9AAB", ["9AAB"]))
    corpus = {"9AAA": corpus_entry("Cryo-EM structure of LY354740-bound mGlu2", [drug, clr, CA]),
              "9AAB": corpus_entry("Cryo-EM structure of apo mGlu2", [clr])}
    obs = {(r["pdb"], r["id"]): r for r in ligands.observe(db, corpus, {"40F": drug, "CLR": clr, "CA": CA})}
    assert obs["9AAA", "40F"]["named"] and obs["9AAA", "40F"]["how"] == "sole"   # the only compound left for the code
    assert obs["9AAA", "40F"]["said"] == "guess" and not obs["9AAA", "40F"]["title"]   # the title's code: only a guess
    assert not obs["9AAA", "CLR"]["named"] and not obs["9AAA", "CA"]["named"]
    assert not obs["9AAB", "CLR"]["named"] and obs["9AAB", "CLR"]["apo"]
    doc = ligands.summarise(list(obs.values()), {"40F": drug})
    assert doc["comps"]["40F"] == {"name": drug["name"], "class": "drug", **tally(1, 1, 0, 0, 1, 1)}
    assert doc["classes"]["lipid"]["all"] == tally(2, 0, 0, 0)


# ---------------------------------------------------------------- classes

def test_classes_come_from_the_component_itself():
    assert ligands.classify(KR0) == "lipid"                                       # 38AQ's fatty acyl fragment
    assert ligands.classify({"name": "3-decanoyloxypropyl decanoate"}) == "lipid"  # by name alone, offline
    assert ligands.classify(CA) == "ion" and ligands.classify({"name": "CALCIUM ION"}) == "ion"
    assert ligands.classify({"name": "SULFATE ION", "formula": "O4 S"}) == "additive"
    assert ligands.classify({"name": "GLYCEROL", "formula": "C3 H8 O3", "smiles": "C(C(CO)O)O"}) == "additive"
    assert ligands.classify({"name": "DODECYL-BETA-D-MALTOSIDE", "type": "D-saccharide"}) == "detergent"
    assert ligands.classify({"name": "2-acetamido-2-deoxy-beta-D-glucopyranose",
                             "type": "D-saccharide, beta linking"}) == "sugar"
    assert ligands.classify({"name": "PROTOPORPHYRIN IX CONTAINING FE", "formula": "C34 H32 Fe N4 O4"}) == "cofactor"
    assert ligands.classify({"name": "IRON/SULFUR CLUSTER", "formula": "Fe4 S4"}) == "cofactor"
    assert ligands.classify(ATP) == "nucleotide"
    assert ligands.classify({"name": "CHOLESTEROL"}) == "lipid"
    assert ligands.classify({"name": "GLUTAMIC ACID", "type": "L-peptide linking"}) == "amino acid"
    assert ligands.classify(NICOTINE) == "drug"
    # 'heptane' inside a drug's systematic name is no alkane
    assert ligands.classify({"name": "N-(2,4,6-trimethylphenyl)bicyclo[2.2.1]heptane-2-carboxamide"}) == "drug"


# ---------------------------------------------------------------- the decision rule

def test_an_incidental_lipid_fragment_is_omitted_with_a_plain_reason(usage):
    d = ligands.decide(KR0, False)
    assert d["decision"] == "omitted" and d["class"] == "lipid"
    assert "lipid or fatty acyl chain" in d["why"] and "does not mention" in d["why"]
    assert "has not seen this one" in d["why"]
    assert ligands.decide(KR0, True)["decision"] == "named"     # unless the PDB title itself names it


def test_a_drug_mpstruc_has_never_seen_is_named(usage):
    d = ligands.decide({"id": "A1XYZ", "name": "examplamide", "weight": 410.0}, False)
    assert d["decision"] == "named" and d["class"] == "drug" and "small-molecule ligand" in d["why"]


def test_an_ion_is_named_only_when_the_title_mentions_it(usage):
    quiet, said = ligands.decide(CA, False), ligands.decide(CA, True)
    assert quiet["decision"] == "omitted" and "32 of the 610" in quiet["why"]
    assert said["decision"] == "named" and "168 of the 190" in said["why"]
    assert ligands.decide({"id": "CL", "name": "CHLORIDE ION", "formula": "Cl"}, False)["decision"] == "omitted"


def test_a_well_known_component_speaks_for_itself_a_rare_one_goes_by_its_class(usage):
    # Cholesterol: 3 of 723 without a title mention; the lipid class does not matter.
    assert ligands.rate(USAGE, "CLR", "lipid", False)["p"] < 0.01
    # A lipid seen once and named once is still pulled most of the way to its class.
    one = {**USAGE, "comps": {"XYZ": {"class": "lipid", **tally(1, 1, 0, 0)}}}
    assert 0.3 < ligands.rate(one, "XYZ", "lipid", False)["p"] < 0.5
    assert ligands.decide({"id": "XYZ", "name": "x"}, False, one)["decision"] == "omitted"


def test_water_and_unknown_atoms_are_never_named(usage):
    for cid in ("HOH", "UNX", "UNL"):
        assert ligands.decide({"id": cid, "name": "x"}, True)["decision"] == "omitted"


def test_tuning_and_scoring_split_by_family():
    assert ligands._fold("9AAA") in ("train", "tune", "test")
    rows = [{"id": "CA", "cls": "ion", "named": False, "said": "none", "apo": False, "family": "F"},
            {"id": "CA", "cls": "ion", "named": True, "said": "title", "apo": False, "family": "F"},
            {"id": "A1XYZ", "cls": "drug", "named": True, "said": "none", "apo": True, "family": "F"}]
    s = ligands.score(rows, USAGE)
    # the apo title silences the drug the title does not mention: one miss
    assert (s["precision"], s["recall"]) == (1.0, 0.5)


# ---------------------------------------------------------------- in a drafted listing

SUSC = entity("TonB-dependent receptor SusC", ["Q8A1G1"], 1003)
SUSD = entity("Starch-binding protein SusD", ["Q8A1G2"], 551)


def test_38aq_the_lipid_fragment_stays_out_of_the_name():
    # With the learned usage as shipped (mpstruc/ligand_usage.json).
    e = entry("38AQ", "Bacteroides SusCD Unliganded - Closed w/ Plug Domain", [SUSC, SUSD],
              [ligand("KR0", "3-decanoyloxypropyl decanoate"), ligand("CA", "CALCIUM ION")])
    mp = build(e, {"Q8A1G1"})
    assert mp["name"] == "TonB-dependent receptor SusC, unliganded, closed state"
    assert "decano" not in mp["name"]
    notes = {n["id"]: n for n in mp["ligandNotes"]}
    assert notes["KR0"]["decision"] == "omitted" and "lipid" in notes["KR0"]["why"]
    assert notes["CA"]["decision"] == "omitted" and notes["CA"]["name"] == "CALCIUM ION"
    assert set(notes["KR0"]) == {"id", "name", "decision", "why"}


def test_notes_explain_named_and_omitted_ligands(usage):
    e = entry("9AAA", "Receptor X in complex with examplamide and calcium", [entity("Receptor X", ["Q1"])],
              [ligand("A1XYZ", "examplamide"), ligand("CA", "CALCIUM ION"), ligand("CLR", "CHOLESTEROL"),
               ligand("HOH", "WATER")])
    mp = build(e, {"Q1"})
    assert mp["name"] == "Receptor X with bound examplamide and Ca<sup>2+</sup>"
    notes = {n["id"]: n["decision"] for n in mp["ligandNotes"]}
    assert notes == {"A1XYZ": "named", "CA": "named", "CLR": "omitted"}          # water is not even listed


def test_an_apo_title_names_no_ligand_it_does_not_mention(usage):
    drug = ligand("A1XYZ", "examplamide")
    e = entry("9AAA", "Receptor X in the apo state", [entity("Receptor X", ["Q1"])], [drug])
    notes: list = []
    assert listing_ligands(e, notes=notes) == []
    assert notes[0]["decision"] == "omitted" and "apo" in notes[0]["why"]
    assert build(e, {"Q1"})["name"] == "Receptor X, apo form"
    # ... but not when the title itself names the ligand
    e = entry("9AAB", "Receptor X bound to examplamide, apo-like conformation", [entity("Receptor X", ["Q1"])], [drug])
    assert [lab for _, lab, _ in listing_ligands(e)] == ["examplamide"]
    # 'unliganded' stays the title's own word, as curators keep it
    e = entry("9AAC", "Unliganded Receptor X", [entity("Receptor X", ["Q1"])], [drug])
    assert build(e, {"Q1"})["name"] == "Receptor X, unliganded"


def test_a_compound_only_the_title_names_goes_to_the_ligand_left_for_it(usage):
    va = ligand("VA", "(Z)-OCTADEC-11-ENYL ACETATE")
    e = entry("9VQP", "The structure of DmOR67d-DmOrco in the cVA-bound state",
              [entity("Odorant receptor 67d", ["Q1"]), entity("Odorant receptor coreceptor", ["Q2"])], [va])
    assert ligands_in_title(e) == {"VA": "compound"}
    assert facets(e, {"Q1", "Q2"})["bound"] == ["cVA"]
    # 'in the presence of' a setting or a protein partner is no compound
    e = entry("9VQQ", "Receptor X in complex with calmodulin in LMNG", [entity("Receptor X", ["Q1"]),
                                                                           entity("Calmodulin-1", ["P0DP23"], 149)], [va])
    assert ligands_in_title(e) == {}
    # 'betaine' in the transporter's own name is still the bound compound of a 'betaine-bound' title
    e = entry("9W99", "Structure of betaine-bound state of the human betaine/GABA transporter 1",
              [entity("Sodium- and chloride-dependent betaine transporter", ["P48065"])],
              [ligand("BET", "TRIMETHYL GLYCINE")])
    assert ligands_in_title(e) == {"BET": "compound"}


def test_an_abbreviation_in_the_title_counts_as_a_mention():
    e = entry("9LI2", "TGR5 in the DCA-bound state", [entity("G-protein coupled bile acid receptor 1", ["Q1"])],
              [ligand("DXC", "(3ALPHA,5BETA,12ALPHA)-3,12-DIHYDROXYCHOLAN-24-OIC ACID", ["DEOXYCHOLIC ACID"])])
    assert ligands_in_title(e) == {"DXC": "abbreviation"}


# ---------------------------------------------------------------- relist

def test_relist_redrafts_listings_and_keeps_every_other_field(tmp_path, monkeypatch):
    doc = {"generated": "2026-11-13T01:52+00:00", "stats": {"released": 2}, "threshold": 0.5,
           "candidates": [
               {"pdb": "38AQ", "title": "Bacteroides SusCD Unliganded", "placement": "parent", "ownerPdb": None,
                "uniprot": ["Q8A1G1"], "membrane": 0.869, "confidence": 0.56, "tmEvidence": [{"uniprot": "Q8A1G1"}],
                "mp": {"name": "SusC with bound 3-decanoyloxypropyl decanoate, closed state", "species": "B. theta",
                       "description": "cryo-EM structure", "curatorNote": "kept"}},
               {"pdb": "38FE", "title": "SusCD + Maltoheptaose", "placement": "related", "ownerPdb": "38AQ",
                "uniprot": ["Q8A1G1"], "membrane": 0.85,
                "mp": {"name": "SusC", "description": "cryo-EM structure", "relatedLine": "old line: 38FE"}},
               {"pdb": "9ZZZ", "title": "no PDB details", "placement": "review", "ownerPdb": None, "uniprot": [],
                "membrane": 0.2, "mp": {"name": "Unchanged", "description": "x-ray structure"}}],
           "inMpstruc": [{"pdb": "1ABC"}]}
    path = tmp_path / "candidates.json"
    path.write_text(json.dumps(doc, separators=(",", ":")))
    seen = []

    def fake_add_listings(props, db=None, summaries=None):
        seen.extend(set(p) for p in props)
        props[0]["mp"] = {"name": "SusC, unliganded, closed state", "description": "cryo-EM structure",
                          "ligandNotes": [{"id": "KR0", "name": "x", "decision": "omitted", "why": "lipid"}]}
        props[1]["mp"] = {"name": "SusC", "description": "cryo-EM structure"}   # no related line any more

    monkeypatch.setattr(automation, "add_listings", fake_add_listings)
    r = automation.relist(path)
    assert seen == [{"pdb", "placement", "ownerPdb", "uniprot"}] * 3            # only what a draft needs
    assert (r["candidates"], r["redrafted"], r["names"], r["descriptions"], r["relatedLines"]) == (3, 2, 1, 0, 1)
    assert [c["pdb"] for c in r["changed"]] == ["38AQ", "38FE"]
    text = path.read_text()
    assert ", " not in text.split('"mp"')[0] and '": ' not in text          # still compact JSON
    new = json.loads(text)
    strip = lambda d: {**d, "candidates": [{k: v for k, v in c.items() if k != "mp"} for c in d["candidates"]]}
    assert strip(new) == strip(doc)                                          # nothing re-scored or re-placed
    a, b, c = (x["mp"] for x in new["candidates"])
    assert a["name"] == "SusC, unliganded, closed state" and a["ligandNotes"][0]["id"] == "KR0"
    assert a["species"] == "B. theta" and a["curatorNote"] == "kept"         # fields the draft does not set stay
    assert "relatedLine" not in b
    assert c == {"name": "Unchanged", "description": "x-ray structure"}      # no PDB details: old draft kept


# ---------------------------------------------------------------- review findings (lipid words, guesses, protein names)

def lig(rec, description=None):
    """A non-polymer entity carrying its chemical component record, as listing.enrich() attaches it."""
    return {"rcsb_nonpolymer_entity": {"pdbx_description": description or rec["name"]},
            "nonpolymer_comp": {"chem_comp": {"id": rec["id"], "name": rec["name"]}, "record": rec,
                                "rcsb_chem_comp_synonyms": [{"name": s[0], "provenance_source": s[1]}
                                                            for s in rec.get("synonyms") or []]}}


REF = "PDB Reference Data"
UBL = {"id": "UBL", "name": "[(2R)-2-oxidanyl-3-phosphonooxy-propyl] (Z)-octadec-9-enoate", "formula": "C21 H41 O7 P",
       "weight": 436.52, "type": "non-polymer", "smiles": "CCCCCCCCC=CCCCCCCCC(=O)OCC(COP(=O)(O)O)O"}
ACOA = {"id": "3IX", "name": "S-[2-[3-[[(2R)-4-[[[(2R,3S,4R,5R)-5-(6-aminopurin-9-yl)-4-oxidanyl-3-phosphonooxy-oxolan-2-yl]"
                              "methoxy-oxidanyl-phosphoryl]oxy-oxidanyl-phosphoryl]oxy-3,3-dimethyl-2-oxidanyl-butanoyl]amino]"
                              "propanoylamino]ethyl] (5Z,8Z,11Z,14Z)-icosa-5,8,11,14-tetraenethioate",
        "formula": "C41 H66 N7 O17 P3 S", "weight": 1053.986, "type": "non-polymer", "synonyms": [["Arachidonyl-CoA", REF]],
        "smiles": "CCCCCC=CCC=CCC=CCC=CCCCC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)n2cnc3c2ncnc3N)O)OP(=O)(O)O)O"}
PT5 = {"id": "PT5", "name": "[(2R)-1-octadecanoyloxy-3-[oxidanyl-[(1R,2R,3S,4R,5R,6S)-2,3,6-tris(oxidanyl)-4,5-diphosphonooxy-"
                             "cyclohexyl]oxy-phosphoryl]oxy-propan-2-yl] (8Z)-icosa-5,8,11,14-tetraenoate",
       "formula": "C47 H85 O19 P3", "weight": 1047.088, "type": "non-polymer",
       "synonyms": [["Phosphatidylinositol 4,5-bisphosphate", REF], ["PtdIns(4,5)P2", REF]],
       "smiles": "CCCCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OC1C(C(C(C(C1O)OP(=O)(O)O)OP(=O)(O)O)O)O)OC(=O)CCCC=CCC=CCC=CCC=CCCCCC"}
CLR = {"id": "CLR", "name": "CHOLESTEROL", "formula": "C27 H46 O", "weight": 386.654, "type": "non-polymer"}
Y01 = {"id": "Y01", "name": "CHOLESTEROL HEMISUCCINATE", "formula": "C31 H50 O4", "weight": 486.726, "type": "non-polymer"}
SCN = {"id": "SCN", "name": "THIOCYANATE ION", "formula": "C N S", "weight": 58.082, "type": "non-polymer"}
PO4 = {"id": "PO4", "name": "PHOSPHATE ION", "formula": "O4 P", "weight": 94.971, "type": "non-polymer"}
ADP = {"id": "ADP", "name": "ADENOSINE-5'-DIPHOSPHATE", "weight": 427.2, "type": "non-polymer"}
MG = {"id": "MG", "name": "MAGNESIUM ION", "formula": "Mg"}
TFX = {"id": "TFX", "name": "2-[4-(dimethylamino)phenyl]-3,6-dimethyl-1,3-benzothiazol-3-ium", "formula": "C17 H19 N2 S",
       "weight": 283.411, "type": "non-polymer", "synonyms": [["Thioflavin T", REF]],
       "smiles": "Cc1ccc2c(c1)sc([n+]2C)c3ccc(cc3)N(C)C"}
A86 = {"id": "A86", "name": "(3S,3'S,5R,5'R,6S,6'R,8'R)-3,5'-dihydroxy-8-oxo-6',7'-didehydro-5,5',6,6',7,8-hexahydro-5,6-epoxy-"
                             "beta,beta-caroten-3'-yl acetate", "synonyms": [["Fucoxanthin", REF]], "drug": "Fucoxanthin"}
CLA = {"id": "CLA", "name": "CHLOROPHYLL A"}


def test_a_lipid_the_title_names_by_a_family_word_is_mentioned():
    # 20ZX: the PDB knows LPA only by a systematic name; with the learned usage as shipped.
    e = entry("20ZX", "Cryo-EM structure of Gq-coupled LPAR5 in complex with LPA",
              [entity("Lysophosphatidic acid receptor 5", ["Q9H1C0"], 372)], [lig(UBL)])
    assert ligands_in_title(e) == {"UBL": "fit"}
    assert build(e, {"Q9H1C0"})["name"].endswith("with bound LPA")
    # 9W63: 'C20:5-CoA bound state'; the cholesterol beside it stays out
    e = entry("9W63", "Cryo-EM structure of C20:5-CoA bound state human ABCD3 in inward-facing conformation",
              [entity("ATP-binding cassette sub-family D member 3", ["P28288"], 914)], [lig(ACOA), lig(CLR)])
    assert ligands_in_title(e) == {"3IX": "fit"}
    name = build(e, {"P28288"})["name"]
    assert "with bound arachidonyl-CoA" in name and "cholesterol" not in name
    # 8JNI: 'in complex with PIP2', not only 'PIP2-bound'
    e = entry("8JNI", "Structure of AE2 in complex with PIP2", [entity("Anion exchange protein 2", ["P04920"], 1241)],
              [lig(PT5), lig({"id": "CL", "name": "CHLORIDE ION", "formula": "Cl"})])
    assert ligands_in_title(e) == {"PT5": "fit"}
    assert build(e, {"P04920"})["name"].endswith("with bound PIP<sub>2</sub>")
    # 'PI(3,5)P2-bound' is read whole; two inositol lipids and one word settle nothing
    e = entry("9Y6P", "Cryo-EM structure of PI(3,5)P2-bound TRPML2", [entity("Mucolipin-2", ["Q1"])], [lig(PT5)])
    assert ligands_in_title(e) == {"PT5": "fit"}
    e["nonpolymer_entities"].append(lig({**PT5, "id": "PIO"}))
    assert ligands_in_title(e) == {}


def test_family_words_go_to_as_many_ligands_as_there_are_words():
    coa = lambda i, name: (i, "lipid", [name])
    two = [coa("A", "3-oxotetracosanoyl-CoA"), coa("B", "Behenoyl-CoA")]
    assert ligands.family_picks(["C22-CoA", "MCoA"], two) == {"A": "C22-CoA", "B": "MCoA"}
    assert ligands.family_picks(["C22-CoA"], two) == {}
    # phosphatidic acid is the phospholipid without a head group; a bare alkyl chain is never the lipid meant
    pa = ("PA", "lipid", ["(2R)-3-(phosphonooxy)propane-1,2-diyl ditetradecanoate"])
    pc = ("PC", "lipid", ["(2S)-3-(hexadecanoyloxy)-2-[(9Z)-octadec-9-enoyloxy]propyl 2-(trimethylammonio)ethyl phosphate"])
    assert ligands.family_picks(["NBD-PA"], [pa, pc, ("C14", "lipid", ["TETRADECANE"])]) == {"PA": "NBD-PA"}
    assert ligands.family_picks(["lipid"], [("C14", "lipid", ["TETRADECANE"])]) == {}
    assert ligands.is_family_word("d16:1 S1P") and ligands.is_family_word("PtdSer") and not ligands.is_family_word("DAMGO")


def test_a_title_compound_is_never_credited_to_a_lipid_by_elimination():
    mor = entity("Mu-type opioid receptor", ["P35372"], 388)
    # 10TM: DAMGO is the five-residue peptide chain, not the cholesterol
    e = entry("10TM", "Cryo-EM structure of MOR bound with DAMGO", [entity("TYR-DAL-GLY-MEA-ETA", [], 5), mor], [lig(CLR)])
    assert ligands_in_title(e) == {}
    assert "cholesterol" not in build(e, {"P35372"})["name"]
    # ... and with no peptide chain either: a lipid is in most structures whatever the title is about
    e = entry("7FIY", "Cryo-EM structure of the tirzepatide-bound human GIPR-Gs complex", [mor], [lig(CLR)])
    assert ligands_in_title(e) == {}
    # a peptide chain takes the title's compound before a small molecule does
    drug = {"id": "A1XYZ", "name": "examplamide", "weight": 410.0}
    e = entry("9AAA", "Cryo-EM structure of MOR bound with DAMGO", [entity("TYR-DAL-GLY-MEA-ETA", [], 5), mor], [lig(drug)])
    assert ligands_in_title(e) == {}
    # 9Z62: 'Na+-bound' is the ion and 'Pi' is phosphate, neither is the cholesterol
    e = entry("9Z62", "Structure of the sodium-dependent phosphate importer SLC34A2, Pi- and Na+-bound",
              [entity("Sodium-dependent phosphate transport protein 2B", ["O95436"])], [lig(CLR), lig(NA), lig(PO4)])
    assert ligands_in_title(e) == {"NA": "name", "PO4": "name"}
    assert [lab for _, lab, _ in listing_ligands(e)] == ["phosphate", "Na<sup>+</sup>"]


def test_an_abbreviation_goes_to_the_one_ligand_it_fits_best():
    # 9VNA: 'with CHS' is cholesterol hemisuccinate (named outright), never cholesterol as well
    e = entry("9VNA", "Cryo-EM structure of TRP melastatin channel in the desensitized state, with CHS",
              [entity("Transient receptor potential cation channel subfamily M member 8", ["Q8R4D5"])],
              [lig(CLR), lig(CA), lig(Y01)])
    assert ligands_in_title(e) == {"Y01": "name"}
    assert build(e, {"Q8R4D5"})["name"].endswith("with bound CHS, desensitized state")
    from mpstruc.listing import _abbreviation_fit
    assert _abbreviation_fit("CHS", ["CHOLESTEROL HEMISUCCINATE"]) > _abbreviation_fit("CHS", ["CHOLESTEROL"]) > 0
    assert _abbreviation_fit("G6P", ["BETA-D-GLUCOSE-6-PHOSPHATE"])      # the digit is in the name
    assert not _abbreviation_fit("CP29", ["CHLOROPHYLL A"])              # a protein's number is not
    # a word of a protein chain's name is no abbreviation of a pigment
    e = entry("9XK8", "PSII-LHCII supercomplex, LHCII trimer with LHCB", [entity("Chlorophyll a-b binding protein LHCB", ["Q1"])],
              [lig({"id": "LUT", "name": "LUTEIN-HYDROXY-CAROTENOID B", "synonyms": [["LUTEIN HYDROXYCAROTENOID B", REF]]})])
    assert ligands_in_title(e) == {}


def test_a_ligand_name_inside_a_protein_name_is_not_a_mention():
    from mpstruc.listing import ligand_said
    said = lambda rec, title: ligand_said(lig(rec), title)
    assert said(ATP, "Porcine ATP synthase with inhibitory protein IF1, DP-state") is None
    assert said(ATP, "ATP-binding cassette transporter in the ATP-bound state") == "ATP"
    assert said(ATP, "Structure of an ATP-sensitive potassium channel") is None
    assert said(ATP, "CNNM4 tetramer with Magnesium and MgATP") == "ATP" and said(ADP, "transporter with MgADPVO4 bound") == "ADP"
    i3p = {"id": "I3P", "name": "D-MYO-INOSITOL-1,4,5-TRIPHOSPHATE"}
    assert said(i3p, "Structure of the IP3 receptor in the closed state") is None
    assert said(i3p, "IP3 receptor bound to IP3 and calcium") == "IP3"
    gly = {"id": "GLY", "name": "GLYCINE"}
    assert said(gly, "Glycine transporter 1 in an occluded state") is None and said(gly, "GlyT1 with glycine") == "glycine"
    # 9PEW: the ATP and ADP of an ATP synthase are not what its title is about (as shipped)
    e = entry("9PEW", "Porcine ATP synthase with inhibitory protein IF1, DP-state",
              [entity("ATP synthase subunit b", ["A0A4X1TA12"], 191)], [lig(ADP), lig(MG), lig(ATP)])
    assert ligands_in_title(e) == {} and listing_ligands(e) == []
    # the same in curators' text
    assert not ligands.mentions("in complex with fucoxanthin chlorophyll a/c-binding protein (FCPI)", A86)
    assert not ligands.mentions("in complex with fucoxanthin chlorophyll-binding protein II (FCPII)", A86)
    assert not ligands.mentions("in complex with alloxanthin chlorophyll a/c-binding protein II (ACPII)", CLA)
    assert ligands.mentions("with bound fucoxanthin", A86) and ligands.mentions("with MgATP bound", ATP)
    assert not ligands.mentions("glycine receptor, open state", gly)


def test_protein_words_in_curator_text_are_no_compounds():
    assert ligands._bound_phrases("with bound LY354740") == ["LY354740"]
    for words in ("With", "5 Lhcr antennae", "auxiliary subunit TARP gamma8"):
        assert ligands._NOT_A_COMPOUND.search(words)
    drug = {"id": "ZEX", "name": "(3R,3'R)-beta,beta-carotene-3,3'-diol", "weight": 568.9}
    db = db_of(mp_entry("5ZGH", "Photosystem I in complex with 5 Lhcr antennae"))
    obs = ligands.observe(db, {"5ZGH": corpus_entry("Cryo-EM structure of the red algal PSI-LHCR", [drug])}, {"ZEX": drug})
    assert not obs[0]["named"]


def test_components_named_x_ion_match_x():
    from mpstruc.listing import _variants, ligand_said
    assert "THIOCYANATE" in _variants("THIOCYANATE ION") and "AZIDE" in _variants("AZIDE ION")
    assert _variants("CALCIUM ION")[-1] == "CALCIUM"       # (ions in listing.IONS never get here)
    # 36EZ, with the usage as shipped
    e = entry("36EZ", "Crystal structure of thiocyanate bound bovine cytochrome c oxidase",
              [entity("Cytochrome c oxidase subunit 1", ["P00396"], 514)], [lig(SCN), lig(MG)])
    assert ligands_in_title(e) == {"SCN": "name"}
    assert build(e, {"P00396"})["name"].endswith("with bound thiocyanate")
    assert ligands.mentions("with bound azide", {"id": "AZI", "name": "AZIDE ION"})
    assert ligands.mentions("hAE2 with bicarbonate", {"id": "BCT", "name": "BICARBONATE ION"})
    af3 = {"id": "AF3", "name": "ALUMINUM FLUORIDE", "formula": "Al F3"}
    assert ligand_said(lig(af3), "NFeoB bound to GDP-AlF3") == "AlF3"


def test_thioflavin_is_a_dye_not_a_flavin_cofactor():
    from mpstruc.listing import title_names
    assert ligands.classify(TFX) == "drug"
    assert ligands.classify({"id": "RBF", "name": "RIBOFLAVIN"}) == "cofactor"
    assert ligands.classify({"id": "FMN", "name": "FLAVIN MONONUCLEOTIDE"}) == "cofactor"
    title = "The cryo-EM structure of ThT bound type1 amyloid beta 42 fibril"
    assert title_names(title) == ["ThT"]
    assert title_names("Cryo-EM structure of Acinetobacter baumannii LolDF bound abaucin (wide cavity)") == []   # protein first
    assert title_names("Structure of GPR3 bound to AF64394") == ["AF64394"]
    e = entry("9XJC", title, [entity("Amyloid-beta precursor protein", ["P05067"], 770)], [lig(TFX)])
    assert build(e, {"P05067"})["name"].endswith("with bound thioflavin T")     # with the usage as shipped


def test_a_guess_is_counted_apart_from_a_mention(usage):
    assert [ligands.situation(h) for h in ("name", "fit", True, "abbreviation", "compound", None, False)] == \
        ["title", "title", "title", "guess", "guess", "none", "none"]
    nuc = {"id": "A1NUC", "name": "adenosine-5'-examplate"}
    assert ligands.decide(nuc, "name")["decision"] == "named"           # nucleotides the title names: 42 of 55
    d = ligands.decide(nuc, "compound")                                 # ... the title only may mean: 6 of 20
    assert d["decision"] == "omitted" and "may be" in d["why"]
    assert ligands.decide({"id": "A1XYZ", "name": "examplamide", "weight": 410.0}, "compound")["decision"] == "named"
    rows = [{"id": "X", "cls": "drug", "named": True, "said": s} for s in ("title", "guess", "guess", "none")]
    assert ligands._tally(rows) == tally(4, 4, 1, 1, 2, 2)
    assert ligands._split(tally(10, 6, 3, 3, 2, 1), "none") == (5, 2)


def test_ligand_notes_are_escaped_text(usage):
    odd = {"id": "A1KB5", "name": "Neu5Ac-alpha2->3Gal-beta1->4Glc & more", "weight": 900.0}
    e = entry("9AAA", "Receptor X", [entity("Receptor X", ["Q1"])], [lig(odd)])
    note = build(e, {"Q1"})["ligandNotes"][0]
    assert note["name"] == "Neu5Ac-alpha2-&gt;3Gal-beta1-&gt;4Glc &amp; more"
    assert "<" not in note["why"] and ">" not in note["why"]


# ---------------------------------------------------------------- the component cache

def test_a_failed_lookup_is_not_cached_as_unknown(tmp_path, monkeypatch):
    from mpstruc import rcsb
    monkeypatch.setattr(rcsb, "CACHE_DIR", tmp_path)
    calls = []
    answer = {"errors": [{"message": "Internal server error"}], "data": None}
    monkeypatch.setattr(rcsb, "_post_json", lambda url, payload: (calls.append(payload["variables"]["ids"]), answer)[1])
    assert ligands.chem_comps(["ATP", "KR0"]) == {}                   # HTTP 200 without data: nothing known ...
    assert not (tmp_path / "chem_comps.json").exists()                # ... and nothing remembered as unknown
    answer = {"data": {"chem_comps": [{"chem_comp": {"id": "ATP", "name": "ADENOSINE-5'-TRIPHOSPHATE",
                                                     "formula_weight": 507.2}}]}}
    got = ligands.chem_comps(["ATP", "KR0"])
    assert calls == [["ATP", "KR0"], ["ATP", "KR0"]]                  # asked again once the PDB answers
    assert got["ATP"]["weight"] == 507.2 and got["KR0"] is None       # the PDB answered: KR0 is not a component
    assert ligands.chem_comps(["ATP", "KR0"]) == got and len(calls) == 2


def test_cached_components_are_used_when_the_pdb_is_down(tmp_path, monkeypatch):
    from mpstruc import listing, rcsb
    monkeypatch.setattr(rcsb, "CACHE_DIR", tmp_path)
    (tmp_path / "chem_comps.json").write_text(json.dumps({"KR0": KR0}))

    def down(url, payload):
        raise OSError("network down")

    monkeypatch.setattr(rcsb, "_post_json", down)
    assert ligands.chem_comps(["KR0", "ZZ9"]) == {"KR0": KR0}
    det = {"38AQ": entry("38AQ", "t", [SUSC], [ligand("KR0", KR0["name"]), ligand("ZZ9", "new thing")])}
    out = listing.enrich(det)
    recs = [n["nonpolymer_comp"].get("record") for n in out["38AQ"]["nonpolymer_entities"]]
    assert recs == [KR0, None]                                        # the cached record is still attached
    assert json.loads((tmp_path / "chem_comps.json").read_text()) == {"KR0": KR0}
