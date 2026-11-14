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


def tally(n, named, n_title, named_title):
    return {"n": n, "named": named, "n_title": n_title, "named_title": named_title}


# A small usage table with the shape of ligand_usage.json.
USAGE = {
    "rule": {"threshold": 0.5, "prior": 2.0},
    "overall": tally(1000, 150, 110, 90),
    "classes": {
        "ion": {"all": tally(500, 50, 60, 45), "rare": tally(40, 8, 8, 8)},
        "lipid": {"all": tally(700, 20, 12, 11), "rare": tally(50, 7, 4, 4)},
        "drug": {"all": tally(300, 220, 120, 114), "rare": tally(200, 150, 70, 68)},
        "nucleotide": {"all": tally(120, 50, 55, 42), "rare": tally(10, 4, 3, 3)},
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
    assert obs["9AAA", "40F"]["title"]                                           # and the title names it the same way
    assert not obs["9AAA", "CLR"]["named"] and not obs["9AAA", "CA"]["named"]
    assert not obs["9AAB", "CLR"]["named"] and obs["9AAB", "CLR"]["apo"]
    doc = ligands.summarise(list(obs.values()), {"40F": drug})
    assert doc["comps"]["40F"] == {"name": drug["name"], "class": "drug", "n": 1, "named": 1, "n_title": 1,
                                   "named_title": 1}
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
    rows = [{"id": "CA", "cls": "ion", "named": False, "title": False, "apo": False, "family": "F"},
            {"id": "CA", "cls": "ion", "named": True, "title": True, "apo": False, "family": "F"},
            {"id": "A1XYZ", "cls": "drug", "named": True, "title": False, "apo": True, "family": "F"}]
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
