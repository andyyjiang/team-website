from mpstruc.automation import (KB, Known, batch_leads, classify, fit_model, in_pool, score,
                                signals, title_tokens)


def _s(core=(), citation=(), pfam=(), sources=(), keywords=(), method="ELECTRON MICROSCOPY", title=""):
    core = set(core)
    return {"uniprot": core, "core": core, "all_uniprot": core, "taxa": {9606}, "pfam": set(pfam),
            "interpro": set(), "opm": set(), "citation": set(citation), "sources": set(sources),
            "keywords": set(keywords), "method": method, "title": title, "n_protein_entities": 1,
            "released": "2020-01-01"}


def _kb():
    kb = KB(entries=[])
    rows = [
        # parent A (P1) with a daughter (same protein, other paper) and a related entry (same paper as A)
        Known("1AAA", "G01-S001-P0001", "parent", "G01-S001", "G01-S001-P0001", "G01-S001-P0001", "2001",
              _s(["P1"], ["doi-a"], ["PF1"])),
        Known("2AAA", "G01-S001-P0001-D001", "daughter", "G01-S001", "G01-S001-P0001", "G01-S001-P0001-D001", "2005",
              _s(["P1"], ["doi-b"], ["PF1"])),
        Known("3AAA", "G01-S001-P0001-D001-R001", "related", "G01-S001", "G01-S001-P0001", "G01-S001-P0001-D001",
              "2005", _s(["P1"], ["doi-b"], ["PF1"])),
        # parent B in another subgroup: a complex of P2 + P3
        Known("1BBB", "G03-S010-P0001", "parent", "G03-S010", "G03-S010-P0001", "G03-S010-P0001", "2010",
              _s(["P2", "P3"], ["doi-c"], ["PF2"])),
    ]
    for k in rows:
        kb.entries.append(k)
        kb.by_pdb.setdefault(k.pdb, []).append(k)
        kb.by_key[k.key] = k
        for u in k.s["core"]:
            kb.by_uniprot.setdefault(u, []).append(k)
        for c in k.s["citation"]:
            kb.by_citation.setdefault(c, []).append(k)
        for f in [f"pfam:{x}" for x in k.s["pfam"]] + [f"up:{u}" for u in k.s["core"]]:
            kb.feat_sg.setdefault(f, __import__("collections").Counter())[k.sg] += 1
        kb.sg_count[k.sg] += 1
    kb.sg_names = {"G01-S001": "Cyclooxygenases", "G03-S010": "Complexes"}
    kb.sg_group = {"G01-S001": "Monotopic", "G03-S010": "Alpha-helical transmembrane"}
    return kb


def test_same_paper_becomes_related_to_where_that_paper_is_listed():
    p = classify(_kb(), "9NEW", _s(["P1"], ["doi-b"]), {}, {})
    assert (p.placement, p.owner_key) == ("related", "G01-S001-P0001-D001")


def test_same_protein_new_paper_becomes_daughter_of_its_parent():
    p = classify(_kb(), "9NEW", _s(["P1"], ["doi-new"]), {}, {})
    assert (p.placement, p.owner_pdb, p.subgroup) == ("daughter", "1AAA", "G01-S001")


def test_fusion_partner_does_not_change_the_protein():
    from mpstruc.automation import summarise
    entry = {"rcsb_primary_citation": {"pdbx_database_id_DOI": "doi-new"}, "polymer_entities": [{
        "entity_poly": {"rcsb_entity_polymer_type": "Protein"},
        "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["P1", "P00720"]},  # T4 lysozyme fusion
        "rcsb_entity_source_organism": [{"ncbi_taxonomy_id": 9606}],
        "rcsb_polymer_entity_annotation": [{"type": "OPM", "annotation_id": "Transmembrane.Alpha.x.y"}],
    }]}
    s = summarise(entry)
    assert s["core"] == {"P1"} and s["uniprot"] == {"P1", "P00720"}
    p = classify(_kb(), "9NEW", s, {}, {})
    assert (p.placement, p.owner_pdb) == ("daughter", "1AAA")


def test_partial_protein_overlap_goes_to_review():
    p = classify(_kb(), "9NEW", _s(["P2"], ["doi-new"]), {}, {})
    assert p.placement == "review" and p.owner_pdb == "1BBB"


def test_new_protein_becomes_parent_in_best_family():
    p = classify(_kb(), "9NEW", _s(["P9"], ["doi-new"], ["PF2"]), {}, {})
    assert (p.placement, p.subgroup) == ("parent", "G03-S010")


def test_no_evidence_goes_to_review():
    p = classify(_kb(), "9NEW", _s(["P9"], ["doi-new"]), {}, {})
    assert p.placement == "review" and p.subgroup is None


def test_same_paper_batch_lists_under_lowest_id_sharing_a_protein():
    cands = {"9ZZB": _s(["P9"], ["doi-x"]), "9ZZA": _s(["P9", "P8"], ["doi-x"]), "9ZZC": _s(["P7"], ["doi-x"])}
    leads = batch_leads(cands)
    p = classify(_kb(), "9ZZB", cands["9ZZB"], leads, {})
    assert (p.placement, p.owner_pdb) == ("related", "9ZZA")
    # A different protein from the same paper is not lumped in.
    assert classify(_kb(), "9ZZC", cands["9ZZC"], leads, {}).placement != "related"
    # The lead itself is not related to anything.
    assert classify(_kb(), "9ZZA", cands["9ZZA"], leads, {}).placement != "related"


def test_signals_and_pool():
    kb = _kb()
    sig = signals(kb, _s(["P1"], ["doi-new"], keywords=["KW-0812"], title="Structure of a sodium channel"))
    assert sig["kw_tm"] == 1 and sig["same_protein"] == 1 and sig["m:em"] == 1 and sig["t:sodium"] == 1
    assert in_pool(sig)
    assert not in_pool(signals(kb, _s(["P9"], ["doi-new"], method="X-RAY DIFFRACTION")))
    assert "t:of" not in title_tokens("Structure of the thing")


def test_model_learns_direction():
    X = [{"a": 1.0}] * 20 + [{"b": 1.0}] * 20
    y = [1] * 18 + [0] * 2 + [0] * 18 + [1] * 2
    m = fit_model(X, y, min_count=1)
    assert score({"a": 1.0}, m) > 0.5 > score({"b": 1.0}, m)
