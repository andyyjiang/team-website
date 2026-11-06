"""Drafted mpstruc listings: names, related lines that name only what differs,
and mpstruc's sub/superscripts (mpstruc/listing.py)."""
from mpstruc.listing import (base_name, build, differences, facets, family_base, ligands, pdb_name,
                             title_facets, title_names, _clean_name, _distinct_title, _esc, _merge, _said)


def entity(desc, uniprot=(), length=400, mutation=None, gene=None):
    return {"entity_poly": {"rcsb_entity_polymer_type": "Protein", "rcsb_sample_sequence_length": length},
            "rcsb_polymer_entity": {"pdbx_description": desc, "pdbx_mutation": mutation},
            "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": list(uniprot)},
            "rcsb_entity_source_organism": [{"ncbi_scientific_name": "Homo sapiens",
                                             "taxonomy_lineage": [{"name": "Eukaryota"}],
                                             "rcsb_gene_name": [{"value": gene}] if gene else None}]}


def ligand(cid, name, synonyms=()):
    return {"nonpolymer_comp": {"chem_comp": {"id": cid, "name": name},
                                "rcsb_chem_comp_synonyms": [{"name": s, "provenance_source": "PDB Reference Data"}
                                                            for s in synonyms]}}


def entry(pdb, title, entities, ligands_=(), res=3.1, method="ELECTRON MICROSCOPY"):
    return {"rcsb_id": pdb, "struct": {"title": title}, "exptl": [{"method": method}],
            "rcsb_entry_info": {"resolution_combined": [res]}, "polymer_entities": list(entities),
            "nonpolymer_entities": list(ligands_), "rcsb_primary_citation": {}}


def complex_i(pdb, title, res=3.1):
    # Two subunits standing in for complex I, with the same cofactor-free ligand set in every map.
    return entry(pdb, title, [entity("NADH-ubiquinone oxidoreductase chain 1", ["P03887"]),
                              entity("NADH dehydrogenase flavoprotein 1", ["P25708"])],
                 [ligand("ADP", "ADENOSINE-5'-DIPHOSPHATE")], res)


CORE_CI = {"P03887", "P25708"}


def test_related_lines_name_only_what_differs_complex_i():
    lead = complex_i("9R12", "Bovine complex I at pH 6.0 in the active state")
    rel = [("9R13", complex_i("9R13", "Bovine complex I at pH 6.0 in the deactive state (consensus)", 3.1)),
           ("9R14", complex_i("9R14", "Bovine complex I at pH 9.2 in the active state", 2.9)),
           ("9R15", complex_i("9R15", "Bovine complex I at pH 9.2 in the deactive state (consensus)", 3.3))]
    m = build(lead, CORE_CI, rel)
    lines = m["description"].split("<br /> ")
    assert lines[0] == "cryo-EM structure"
    # ADP is bound in every map, so no line mentions it; pH and state differ, so every line gives them.
    assert lines[1] == "pH 6.0, deactive state, consensus, 3.10 &Aring;: 9R13"
    assert lines[2] == "pH 9.2, active state, 2.90 &Aring;: 9R14"
    assert lines[3] == "pH 9.2, deactive state, consensus, 3.30 &Aring;: 9R15"
    assert "ADP" not in m["description"]


RECEPTOR = entity("5-hydroxytryptamine receptor 2A", ["P28223"], gene="HTR2A")
G_ALPHA = entity("Guanine nucleotide-binding protein G(s) subunit alpha isoforms short", [], 380)
G_BETA = entity("Guanine nucleotide-binding protein G(I)/G(S)/G(T) subunit beta-1", ["P62873"], 340)
NB35 = entity("Nanobody 35", [], 130)


def gpcr(pdb, title, lig):
    return entry(pdb, title, [RECEPTOR, G_ALPHA, G_BETA, NB35], [lig])


def test_shared_partners_left_out_only_the_ligand_differs():
    lead = gpcr("9AAA", "5-HT2AR bound to LSD in complex with Gs", ligand("7LD", "(8alpha)-N,N-diethyl-6-methyl"))
    rel = [("9AAB", gpcr("9AAB", "5-HT2AR bound to IHCH-1904 in complex with Gs", ligand("A1EYI", "(2S)-2-[(4,5)-x]"))),
           ("9AAC", gpcr("9AAC", "5-HT2AR bound to psilocin in complex with Gs", ligand("PSI", "PSILOCIN")))]
    d = differences([("9AAA", lead)] + rel, {"P28223"})
    # The G protein and the nanobody are in every structure: just 'with bound X'.
    assert d["9AAB"] == "with bound IHCH-1904"  # a new compound takes the name the title gives it
    assert d["9AAC"] == "with bound psilocin"


def test_lead_condition_goes_in_the_name_with_the_family_name_and_markup():
    lead = gpcr("9AAA", "5-HT2AR bound to LSD in complex with mini-Gq and scFv16",
                ligand("7LD", "(8alpha)-N,N-diethyl-6-methyl", ["Lysergic acid diethylamide"]))
    fam = [("5-HT<sub>2A</sub> serotonin receptor with bound risperidone", "2018"),
           ("5-HT<sub>2A</sub> serotonin receptor - G<sub>i</sub> complex", "2025")]
    m = build(lead, {"P28223"}, (), {"family_names": fam})
    assert m["name"] == "5-HT<sub>2A</sub> serotonin receptor - miniG<sub>q</sub> complex with bound LSD"
    assert m["description"] == "cryo-EM structure"


def test_pdb_name_gets_mpstruc_markup_and_gene_symbol():
    e = entry("9AAD", "Structure of apo MgtE", [entity("Mg2+ transporter MgtE", ["Q1"], gene="mgtE")])
    m = build(e, {"Q1"})
    assert m["name"] == "Mg<sup>2+</sup> transporter MgtE, apo form"
    e = entry("9AAE", "Structure of MCT2", [entity("Monocarboxylate transporter 2", ["Q2"], gene="SLC16A7")])
    assert pdb_name(e["polymer_entities"]) == "SLC16A7 Monocarboxylate transporter 2"


def test_bibliography_title_stays_as_published():
    e = entry("9AAF", "x", [entity("Channel", ["Q1"])])
    e["rcsb_primary_citation"] = {"title": "Gating of the K+ channel by Ca2+", "rcsb_authors": ["Li, X."]}
    # mpstruc adds no markup to paper titles (none of 803 titles for 2024-2026 papers has any).
    assert build(e, {"Q1"})["bibliography"]["title"] == "Gating of the K+ channel by Ca2+"


def test_nothing_differs_falls_back_to_title_words_then_additional_structure():
    a = entry("9BBA", "Structure of transporter X", [entity("Transporter X", ["Q1"])])
    b = entry("9BBB", "Structure of transporter X, alternative conformation", [entity("Transporter X", ["Q1"])])
    c = entry("9BBC", "Structure of transporter X", [entity("Transporter X", ["Q1"])])
    d = differences([("9BBA", a), ("9BBB", b), ("9BBC", c)], {"Q1"})
    assert d["9BBB"] == "alternative conformation"
    assert d["9BBC"] == "additional structure"


def test_mutants_partners_and_peptides():
    wt = entry("9CCA", "Human PIEZO1-MDFIC", [entity("Piezo-type mechanosensitive ion channel component 1", ["Q92508"],
                                                   gene="PIEZO1"),
                                            entity("MyoD family inhibitor domain-containing protein", ["Q9P1T7"],
                                                   gene="MDFIC")])
    mut = entry("9CCB", "Human PIEZO1-A1988V", [entity("Piezo-type mechanosensitive ion channel component 1",
                                                      ["Q92508"], mutation="A1988V", gene="PIEZO1")])
    d = differences([("9CCA", wt), ("9CCB", mut)], {"Q92508"})
    assert d["9CCB"] == "A1988V mutant"
    tox = entry("9CCC", "Human Slo1-Charybdotoxin complex", [
        entity("Calcium-activated potassium channel subunit alpha-1", ["Q12791"], gene="SLO"),
        entity("Potassium channel toxin alpha-KTx 1.1", ["P13487"], 37)], [ligand("K", "POTASSIUM ION")])
    f = facets(tox, {"Q12791"})
    # The toxin is a short chain: 'with bound', named as the title names it (in mpstruc's lower case);
    # 'potassium channel' names the protein, so K+ is not listed as a bound ion.
    assert f["bound"] == ["charybdotoxin"]


def test_title_names_skip_the_protein_itself():
    # (Not one of the 19 entries curators added 2-15 Oct 2026, which are kept for checking the drafts.)
    assert title_names("Cryo-EM structure of mGlu2 bound to agonist LY354740 and PAM JNJ-46281222",
                       {"mGlu2"}) == ["LY354740", "JNJ-46281222"]
    assert title_names("Human Slo1-Charybdotoxin complex", {"SLO"}) == ["Charybdotoxin"]
    assert title_names("structure of MCT2-embigin-AR-C155858 complex", {"MCT2"}) == ["embigin", "AR-C155858"]


def test_ligand_filters_and_spellings():
    e = entry("9DDA", "Cdr1 with AMPPNP", [entity("Cdr1", ["Q1"])],
              [ligand("ANP", "PHOSPHOAMINOPHOSPHONIC ACID-ADENYLATE ESTER"), ligand("MG", "MAGNESIUM ION"),
               ligand("LMT", "DODECYL-BETA-D-MALTOSIDE"), ligand("CLA", "CHLOROPHYLL A"),
               ligand("PAX", "Paxilline")])
    labels = [lab for _, lab, _ in ligands(e)]
    # The title's spelling of AMP-PNP; ions, detergents and cofactors the title does not name are left out;
    # drug names in mpstruc's lower case.
    assert labels == ["AMPPNP", "paxilline"]


def test_base_name_and_family_base():
    assert base_name("C5a anaphylatoxin chemotactic receptor 1 (C5aR) - G<sub>i</sub> complex, with bound C5a") == \
        "C5a anaphylatoxin chemotactic receptor 1 (C5aR)"
    assert base_name("Monocarboxylate transporter 2 (MCT2, SLC16 family)") == "Monocarboxylate transporter 2 (MCT2, SLC16 family)"
    assert base_name("InsP<sub>3</sub>R3 Inositol-1,4,5-trisphosphate receptor, apo state") == \
        "InsP<sub>3</sub>R3 Inositol-1,4,5-trisphosphate receptor"
    assert base_name("Metabotropic Glutamate Receptor 8 (mGlu8), apo form") == "Metabotropic Glutamate Receptor 8 (mGlu8)"
    assert base_name("Electron Transport Chain (ETC) Super Complex III<sub>2</sub>IV") == \
        "Electron Transport Chain (ETC) Super Complex III<sub>2</sub>IV"
    fam = [("glucagon receptor (GCGR)", "2013"), ("glucagon receptor (GCGR)- G<sub>s</sub> complex, apo form", "2024"),
           ("Glucagon receptor (GCGR)-G<sub>s</sub> complex with bound Peptide 15", "2025")]
    assert family_base(fam) == "Glucagon receptor (GCGR)"


def test_related_line_against_the_entry_it_goes_under():
    owner = entry("8AAA", "Transporter T, apo state", [entity("Transporter T", ["Q1"])])
    known = entry("8AAB", "Transporter T bound to substrate S1", [entity("Transporter T", ["Q1"])],
                  [ligand("SB1", "Substrate S1")], 2.8)
    new = entry("9EEA", "Transporter T bound to inhibitor ZX-123", [entity("Transporter T", ["Q1"])],
                [ligand("A1ZX", "(2R)-2-[4-(1,2)-y]")], 2.95)
    m = build(new, {"Q1"}, (), {"group": [("8AAA", owner), ("8AAB", known)]})
    assert m["relatedLine"] == "with bound ZX-123, 2.95 &Aring;: 9EEA"


def test_enrich_merges_extra_fields_by_entity_name():
    det = entry("9FFA", "t", [entity("Channel", ["Q1"])], [ligand("PAX", "PAXILLINE")])
    det["polymer_entities"][0]["rcsb_polymer_entity"]["pdbx_mutation"] = None
    extra = {"polymer_entities": [{"rcsb_polymer_entity": {"pdbx_description": "Channel", "pdbx_mutation": "E71A"},
                                   "entity_poly": {"rcsb_sample_sequence_length": 99},
                                   "rcsb_entity_source_organism": [{"rcsb_gene_name": [{"value": "KCSA"}]}]}],
             "nonpolymer_entities": [{"rcsb_nonpolymer_entity": {"pdbx_description": "Paxilline"},
                                      "nonpolymer_comp": {"chem_comp": {"id": "PAX"},
                                                          "rcsb_chem_comp_synonyms": [{"name": "Pax", "provenance_source": "Author"}]}}]}
    _merge(det, extra)
    e = det["polymer_entities"][0]
    assert e["rcsb_polymer_entity"]["pdbx_mutation"] == "E71A" and e["entity_poly"]["rcsb_sample_sequence_length"] == 99
    assert e["rcsb_entity_source_organism"][0]["rcsb_gene_name"] == [{"value": "KCSA"}]
    assert det["nonpolymer_entities"][0]["rcsb_nonpolymer_entity"]["pdbx_description"] == "Paxilline"


def test_structure_that_only_lacks_a_ligand_says_what_it_has():
    both = entry("9GGA", "MOR bound with DAMGO and BMS986187", [entity("Mu-type opioid receptor", ["P35372"])],
                 [ligand("DAM", "DAMGO"), ligand("A1BM", "(3R)-x-[1,2]-y")])
    one = entry("9GGB", "MOR bound with DAMGO", [entity("Mu-type opioid receptor", ["P35372"])], [ligand("DAM", "DAMGO")])
    assert differences([("9GGA", both), ("9GGB", one)], {"P35372"})["9GGB"] == "with bound DAMGO"



# ---------------------------------------------------------------- review fixes

def channel(pdb, title, mutation=None, res=3.0, ligands_=()):
    return entry(pdb, title, [entity("Transient receptor potential cation channel subfamily V member 3", ["Q8NET8"],
                                     mutation=mutation, gene="TRPV3")], ligands_, res)


def test_siblings_never_share_a_line_when_their_titles_differ():
    lead = channel("9PJ6", "Inactivated-state structure of wild-type human TRPV3 purified in GDN detergent")
    rel = [("13LQ", channel("13LQ", "Structure of human TRPV3-G568D mutant in the 4-fold symmetrical open state", "G568D")),
           ("13LT", channel("13LT", "Structure of human TRPV3-G568D mutant in the 2-fold symmetrical open state", "G568D")),
           ("13LN", channel("13LN", "Structure of human TRPV3-Q580P mutant in the inactivated state", "Q580P")),
           ("13LX", channel("13LX", "Structure of human TRPV3-G568D mutant, open state, dataset alpha", "G568D")),
           ("13LY", channel("13LY", "Structure of human TRPV3-G568D mutant, open state, dataset beta", "G568D"))]
    d = differences([("9PJ6", lead)] + rel, {"Q8NET8"})
    assert d["13LQ"] == "G568D mutant, 4-fold symmetrical open state"
    assert d["13LT"] == "G568D mutant, 2-fold symmetrical open state"
    assert d["13LN"] == "Q580P mutant, inactivated state"  # curators' word, no longer lost
    # Same mutation and state: the title words that differ are added.
    assert d["13LX"].endswith("alpha") and d["13LY"].endswith("beta")
    titles = {p: det["struct"]["title"] for p, det in rel}
    lines = [d[p] for p in titles]
    assert len(set(lines)) == len(lines)


def test_state_words_numbered_conformations_and_time_points():
    assert title_facets("hOCT1 in inward facing partially open 2 conformation")["state"] == \
        ["inward facing partially open state 2"]
    assert title_facets("Primed-state naloxone-mu opioid receptor-Gi complex")["state"] == ["primed state"]
    for word in ("latent", "unlatched", "engaged", "inactivated"):
        assert title_facets(f"Receptor X in the {word} state")["state"] == [f"{word} state"]
    assert title_facets("Tetrameric 5-HT3A receptor in Salipro (apo, asymmetric)")["state"][-1] == "asymmetric"
    assert title_facets("Na+-PPase 300-seconds post reaction initiation")["delay"] == ["300 seconds"]
    assert title_facets("Na+-PPase 0-60-seconds post reaction initiation")["delay"] == ["0-60 seconds"]
    assert "state" not in title_facets("FepE closed 9 subunit complex C9 symmetry") or \
        title_facets("FepE closed 9 subunit complex")["state"] == ["closed state"]


def test_title_word_fallback_drops_measurements_and_filler():
    t = ["Photosystem II dark-stable (0F, S1 rich) at 1.97 Angstrom obtained from conventional data at 9.5 keV",
         "Photosystem II one-flash (1F, S2 rich) at 2.10 Angstrom obtained from anomalous data at 9.5 keV",
         "Photosystem II dark-stable (0F, S1 rich) at 2.16 Angstrom obtained from anomalous data at 7.1 keV"]
    assert _distinct_title(t, 0, "") == "dark-stable 0F S1, conventional"
    assert _distinct_title(t, 1, "") == "one-flash 1F S2, anomalous"
    assert _distinct_title(["Flagellar motor in situ, inner ring", "Flagellar motor in situ, outer ring"], 0, "") == "inner"
    assert _distinct_title(["KCNQ2/3 with 2:2 stoichiometry", "KCNQ2/3 with 3:1 stoichiometry"], 0, "") == "2:2"
    assert _distinct_title(["ABCB4 in presence of elacridar", "ABCB4 in presence of zosuquidar"], 0, "") == "elacridar"
    assert _distinct_title(["BR light dataset at 2493 GW/cm2", "BR light dataset at 342 GW/cm2"], 0, "") == "2493 GW/cm2"


def test_partner_labels_are_never_plain_words():
    pep = {"rcsb_polymer_entity": {"pdbx_description": "N-terminal fragment of the engineered fractalkine CX3CL1.44"}}
    assert _said(pep, "Structure of the X-Y complex in the E-state") is None
    nup = {"rcsb_polymer_entity": {"pdbx_description": "Nuclear pore complex protein Nup155"}}
    assert _said(nup, "Integrative structure of the human nuclear pore complex nuclear basket") is None
    assert _said(nup, "Structure of the Nup155 complex") == "Nup155"
    # A 9XH2-like lead: the G protein named from the title (Gq) and from the chains (Gi) of a chimera is one
    # partner, and the peptide keeps a real name, so related lines give only what differs.
    rec = entity("US28 GPCR", ["P69332"], 350)
    g_a = entity("Guanine nucleotide-binding protein G(i) subunit alpha-1", [], 354)
    lig = entity("N-terminal fragment of the engineered fractalkine CX3CL1.44", ["P78423"], 80)
    lead = entry("9XH2", "Structure of US28 bound to CX3CL1 in complex with Gq", [rec, g_a, lig])
    rel = [("9XIX", entry("9XIX", "US28-CX3CL1 complex in the C-state", [rec, g_a, lig])),
           ("9XIY", entry("9XIY", "US28-CX3CL1 complex in the T2C-state with GDP", [rec, g_a, lig],
                          [ligand("GDP", "GUANOSINE-5'-DIPHOSPHATE")]))]
    m = build(lead, {"P69332"}, rel)
    assert "the" not in m["name"].split() and "G<sub>i</sub>" not in m["description"]
    assert "CX3CL1" in m["name"]


def test_names_from_pdb_entities_are_clean():
    assert _clean_name("2-oxoglutarate receptor 1,2-oxoglutarate receptor 1,2-oxoglutarate receptor 1,Fusion protein,"
                       "Maltose/maltodextrin-binding periplasmic protein") == "2-oxoglutarate receptor 1"
    assert _clean_name("Beta-2 adrenergic receptor,D(1A) dopamine receptor,Engineered oplophorus-luciferin "
                       "2-monooxygenase catalytic subunit") == "Beta-2 adrenergic receptor, D(1A) dopamine receptor"
    assert _clean_name("Inositol 1,4,5-trisphosphate receptor type 1") == "Inositol 1,4,5-trisphosphate receptor type 1"
    # PDB charges become superscripts; no gene symbol before a list of fused parts.
    e = entry("9AAG", "Structure of apo CLC-3", [entity("H(+)/Cl(-) exchange transporter 3", ["Q1"], gene="CLCN3")])
    assert build(e, {"Q1"})["name"] == "CLCN3 H<sup>+</sup>/Cl<sup>-</sup> exchange transporter 3, apo form"
    e = entry("9AAH", "x", [entity("Beta-2 adrenergic receptor,D(1A) dopamine receptor", ["Q1"], gene="ADRB2")])
    assert pdb_name(e["polymer_entities"]) == "Beta-2 adrenergic receptor, D(1A) dopamine receptor"
    # Two state words become one 'X, Y state'; the same position is never mutated twice.
    e = entry("9AAI", "GlyT2 in the inward-open, substrate-free state", [entity("Glycine transporter 2", ["Q1"])])
    assert build(e, {"Q1"})["name"] == "Glycine transporter 2, inward-open, substrate-free state"
    e = entry("9AAJ", "RyR1 I4657M/G4819E mutant", [entity("Ryanodine receptor 1", ["Q1"], mutation="C4657M,G4819E")])
    assert facets(e, {"Q1"})["mut"] == ["I4657M", "G4819E"]


def test_ligand_names_as_curators_write_them():
    e = entry("9DDB", "TIR1 dimer with auxin", [entity("TIR1", ["Q1"])],
              [ligand("IAC", "1H-INDOL-3-YLACETIC ACID", ["INDOLE ACETIC ACID"])])
    # The title's own word wins over the chemical name.
    assert facets(e, {"Q1"})["bound"] == ["auxin"]
    e = entry("9DDC", "Pump in complex with imidodiphosphate", [entity("Pump", ["Q1"])],
              [ligand("2PN", "IMIDODIPHOSPHORIC ACID")])
    assert [lab for _, lab, _ in ligands(e)] == ["imidodiphosphate"]
    e = entry("9DDD", "Receptor X structure", [entity("Receptor X", ["Q1"])],
              [ligand("NCT", "(S)-(-)-NICOTINE"), ligand("6BA", "N-BENZYL-9H-PURIN-6-AMINE"),
               ligand("FES", "FE2/S2 (INORGANIC) CLUSTER")])
    labels = [lab for _, lab, _ in ligands(e)]
    assert "nicotine" in labels and "N-benzyl-9H-purin-6-amine" in labels  # locants keep their capitals
    e = entry("9DDE", "Photosystem II", [entity("Photosystem II protein D1", ["Q1"])],
              [ligand("LNL", "alpha-linolenic acid"), ligand("ZZM", "(2S)-x-propyl octadecanoate",
                                                          ["monogalactosyl-diacylglycerol"])])
    assert ligands(e) == []  # lipids and fatty acids, by any of their names


def test_peptides_by_sequence_and_large_complex_subunits():
    pep = entity("LYS-LYS-PRO-TYR-ILE-LEU", [], 6)
    rec = entity("Neurotensin receptor type 1", ["P30989"], gene="NTSR1")
    e = entry("9VAA", "NTSR1-Gi complex", [rec, pep])
    assert facets(e, {"P30989"})["bound"] == ["peptide"]
    e = entry("9VAB", "NTSR1 bound to neurotensin", [rec, pep])
    assert facets(e, {"P30989"})["bound"] == ["neurotensin"]  # the title's name for it
    psii = [entity(f"Photosystem II reaction center protein {x}", [f"P{i}"], 300) for i, x in enumerate("ABCD")]
    e = entry("9XK3", "Photosystem II supercomplex", psii + [entity("Os02g0581100 protein", ["Q9"], 90)])
    f = facets(e, {"P0", "P1", "P2", "P3"})
    assert f["bound"] == [] and f["cplx"] == []
    # Another subunit of the same ABC transporter is not 'in complex with'.
    abc = entry("9Q2N", "PmtCD in nanodisc", [entity("Phenol-soluble modulin export ABC transporter permease subunit PmtD",
                                                     ["Q2FWX0"]), entity("ABC transporter ATP-binding protein", ["Q2FWW9"])])
    assert facets(abc, {"Q2FWX0"})["cplx"] == []


def test_g_protein_subtype_spelling_is_not_a_difference():
    def mor(pdb, title, g_desc):
        return entry(pdb, title, [entity("Mu-type opioid receptor", ["P35372"]), entity(g_desc, [], 354)])
    lead = mor("9PXV", "Nucleotide-free naloxone-mu opioid receptor Gi1 complex",
               "Guanine nucleotide-binding protein G(i) subunit alpha-1")
    rel = [(p, mor(p, f"{s}-state naloxone-mu opioid receptor-Gi complex",
                   "Guanine nucleotide-binding protein G(i) subunit alpha-1"))
           for p, s in (("9PXW", "Latent"), ("9PXX", "Unlatched"), ("9PXY", "Primed"))]
    d = differences([("9PXV", lead)] + rel, {"P35372"})
    assert d == {"9PXV": "nucleotide-free state", "9PXW": "latent state", "9PXX": "unlatched state",
                 "9PXY": "primed state"}


def test_markup_from_pdb_text_is_balanced_and_escaped():
    assert _esc("Ca<sup>2+</sup> and x<sub>2</sub>") == "Ca<sup>2+</sup> and x<sub>2</sub>"
    assert _esc("foo<sub>bar") == "foobar"  # an unclosed tag is dropped, never restored
    assert _esc("a < b > c & d") == "a &lt; b &gt; c &amp; d"
    bad = "<img src=x onerror=alert(1)>"
    lead = entry("9AAA", f"Channel X <sub>bound to foo", [entity(f"Channel X {bad}", ["P1"])],
                 [ligand("ZZZ", "foo<sub>bar"), ligand("ZZY", "a<b & c>d")])
    rel = entry("9AAB", f"Channel X with </sup>{bad} & Ca2+", [entity("Channel X", ["P1"])],
                [ligand("QQQ", "x</sub><img src=y onerror=alert(2)>")])
    m = build(lead, {"P1"}, [("9AAB", rel)], {"group": [("9AAB", rel)]})
    for text in (m["name"], m["description"], m["relatedLine"]):
        assert "<img" not in text and "<script" not in text
        assert text.count("<sub>") == text.count("</sub>") and text.count("<sup>") == text.count("</sup>")
