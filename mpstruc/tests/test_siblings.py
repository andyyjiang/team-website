"""Same-paper structures released outside the weekly window (mpstruc/siblings.py)."""
import json

from mpstruc import automation, rcsb, siblings
from mpstruc.automation import KB, batch_leads, classify, mark_siblings

DOI = "10.1038/s41467-026-71212-6"
KEY = f"doi:{DOI}"
PAPER = {"journal_abbrev": "Nat Commun", "pdbx_database_id_DOI": DOI, "pdbx_database_id_PubMed": 41888148}
WINDOW = ["2026-05-15", "2026-11-13"]


def _s(core, citation, released):
    core = set(core)
    return {"uniprot": core, "core": core, "all_uniprot": core, "taxa": {9606}, "pfam": set(), "interpro": set(),
            "opm": set(), "citation": set(citation), "sources": set(), "keywords": set(), "method": "ELECTRON MICROSCOPY",
            "title": "", "n_protein_entities": 1, "released": released}


# ---------------------------------------------------------------- papers

def test_paper_key_prefers_the_doi_and_falls_back_to_pubmed():
    assert siblings.paper_key("10.1038/S41467-026-71212-6", 41888148) == KEY
    assert siblings.paper_key(None, 41888148) == "pmid:41888148"
    assert siblings.paper_key("", "-1") is None and siblings.paper_key(None, None) is None


def test_papers_come_only_from_published_candidates_of_the_window():
    cands = [{"pdb": "9WFO", "doi": DOI, "pmid": "41888148"},
             {"pdb": "9AAA", "doi": "10.1101/2026.01.01.123456", "awaiting": {"reason": "preprint only"}},
             {"pdb": "9WF2", "doi": DOI, "sibling": {"of": DOI}},  # found by an earlier call: brings nothing in
             {"pdb": "9BBB", "doi": None, "pmid": "123"}, {"pdb": "9CCC", "doi": None, "pmid": None}]
    assert siblings.papers(cands) == {KEY: ["9WFO"], "pmid:123": ["9BBB"]}


def test_one_search_asks_for_dois_and_pubmed_ids_together():
    q = siblings.search_query([KEY, "pmid:123"])
    nodes = q["query"]["nodes"]
    assert q["query"]["logical_operator"] == "or" and q["return_type"] == "entry"
    assert nodes[0]["parameters"] == {"attribute": siblings.DOI_ATTR, "operator": "in", "value": [DOI]}
    assert nodes[1]["parameters"] == {"attribute": siblings.PMID_ATTR, "operator": "in", "value": [123]}
    assert siblings.search_query([KEY])["query"]["type"] == "terminal"


# ---------------------------------------------------------------- discovery

def test_find_entries_batches_matches_by_citation_and_caches(tmp_path):
    cache = tmp_path / "paper_entries.json"
    keys = [KEY, "pmid:123"] + [f"doi:10.1/x{i}" for i in range(siblings.BATCH)]
    calls = []

    def search(q):
        calls.append(q)
        return ["9WF2", "9WFO", "9BBB", "9ODD"] if len(calls) == 2 else []  # (the papers are asked in sorted order)

    cites = {"9WF2": dict(PAPER, pdbx_database_id_DOI=DOI.upper()), "9WFO": PAPER,
             "9BBB": {"pdbx_database_id_DOI": "10.9/other", "pdbx_database_id_PubMed": 123},
             "9ODD": {"pdbx_database_id_DOI": "10.9/unasked"}}
    found, got = siblings.find_entries(keys, search, lambda ids: {p: cites[p] for p in ids}, cache, now=1000.0, pause=0)
    assert len(calls) == 2  # 52 papers, 50 per request
    assert found[KEY] == ["9WF2", "9WFO"] and found["pmid:123"] == ["9BBB"] and found["doi:10.1/x3"] == []
    assert set(got) == {"9WF2", "9WFO", "9BBB", "9ODD"}
    # Within the cache's life nothing is asked again; after it, everything is.
    again, got = siblings.find_entries(keys, lambda q: 1 / 0, None, cache, now=1000.0 + siblings.MAX_AGE_S - 1, pause=0)
    assert again == found and got == {}
    calls.clear()
    siblings.find_entries([KEY], search, lambda ids: {p: cites[p] for p in ids}, cache, now=1000.0 + siblings.MAX_AGE_S + 1, pause=0)
    assert len(calls) == 1


def test_find_entries_survives_a_failed_lookup_and_asks_again_next_time(tmp_path, capsys):
    cache = tmp_path / "paper_entries.json"

    def down(q):
        raise OSError("search is down")

    found, got = siblings.find_entries([KEY], down, lambda ids: {}, cache, now=1.0, pause=0)
    assert found == {} and got == {} and not cache.exists()
    assert "could not look up 1 papers" in capsys.readouterr().out
    found, _ = siblings.find_entries([KEY], lambda q: ["9WF2"], lambda ids: {"9WF2": PAPER}, cache, now=2.0, pause=0)
    assert found == {KEY: ["9WF2"]}


def test_select_takes_what_the_run_and_mpstruc_do_not_have():
    found = {KEY: ["9WF2", "9WF3", "9WF9", "9WFC", "9WFL", "9WFN", "9WFO", "9NEW", "9XXX"]}
    released = {"9WF2": "2026-04-22", "9WF3": "2026-04-01", "9WF9": "2026-04-01", "9WFC": "2026-04-01",
                "9WFL": "2026-04-01", "9WFN": "2026-04-01", "9NEW": "2026-05-15", "9XXX": ""}
    sibs = siblings.select({KEY: ["9WFO"]}, found, released, WINDOW, have={"9WFO"}, listed={"9WF9"},
                           decided={"9WFC"}, unpublished=lambda p: p == "9WFL")
    # 9WFO is the candidate itself, 9WF9 is in mpstruc, 9WFC was decided by a curator, 9WFL's own
    # record has no paper yet, 9NEW was released inside the window, 9XXX has no release date.
    assert sibs == {p: {"paper": KEY, "with": ["9WFO"]} for p in ("9WF2", "9WF3", "9WFN")}


def test_stale_cached_records_are_those_naming_another_paper():
    cached = {"9WF2": {"rcsb_primary_citation": {"pdbx_database_id_DOI": None}}, "9WF3": {"rcsb_primary_citation": PAPER},
              "9WF9": None}
    assert siblings.stale(cached, {"9WF2": PAPER, "9WF3": PAPER, "9WF9": PAPER, "9ZZZ": PAPER}) == ["9WF2", "9WF9"]


def test_discover_drops_stale_cache_entries_and_reads_release_dates(monkeypatch):
    forgot = []
    monkeypatch.setattr(siblings.awaiting, "forget", lambda ids: forgot.append(list(ids)))
    monkeypatch.setattr(siblings.awaiting, "fresh_citations", lambda ids: {p: PAPER for p in ids})
    dates = {"9WF2": "2026-04-22", "9WF3": "2026-04-01", "9OLD": "2026-11-01"}
    monkeypatch.setattr(rcsb, "entry_entities", lambda ids: {p: {
        "rcsb_accession_info": {"initial_release_date": dates[p] + "T00:00:00Z"},
        "rcsb_primary_citation": None if p == "9WF2" else PAPER} for p in ids})
    find = lambda keys: ({KEY: ["9WF2", "9WF3", "9WFO", "9OLD"]}, {"9WF2": PAPER})  # 9WF3 and 9OLD came from the cache
    sibs, counts = siblings.discover([{"pdb": "9WFO", "doi": DOI}], WINDOW, {"9WFO"}, set(), set(), find)
    assert forgot == [["9WF2"]]
    assert sibs == {"9WF2": {"paper": KEY, "with": ["9WFO"], "released": "2026-04-22"},
                    "9WF3": {"paper": KEY, "with": ["9WFO"], "released": "2026-04-01"}}
    assert counts["papers"] == 1 and counts["found"] == 2 and counts["papersWithSiblings"] == 1


# ---------------------------------------------------------------- keeping rule

def test_keep_follows_the_floor_except_for_the_protein_of_a_kept_candidate_of_the_paper():
    sibs = {p: {"paper": KEY, "with": ["9WFO", "9LOW"]} for p in ("9HIG", "9SAM", "9DIF", "9LOS", "9NOS")}
    scores = {"9HIG": 0.7, "9SAM": 0.04, "9DIF": 0.04, "9LOS": 0.04}  # 9NOS could not be scored
    cores = {"9WFO": {"P38435"}, "9LOW": {"Q1"}, "9HIG": {"P9"}, "9SAM": {"P38435"}, "9DIF": {"P38435", "P2"},
             "9LOS": {"Q1"}}
    got = siblings.keep(sibs, scores, 0.1, cores, kept={"9WFO"})  # 9LOW was itself below the floor
    assert got == {"9HIG": "score", "9SAM": "same protein as 9WFO"}


# ---------------------------------------------------------------- lead and notes

def test_the_lead_of_a_paper_with_siblings_is_its_lowest_id_whatever_the_release_dates():
    kb = KB(entries=[])
    cands = {"9WFO": _s(["P38435"], [KEY], "2026-07-15"), "9WF2": _s(["P38435"], [KEY], "2026-04-22"),
             "9WF3": _s(["P38435"], [KEY], "2026-04-01"), "9WFL": _s(["P38435", "O14668"], [KEY], "2026-04-01")}
    leads = batch_leads(cands)
    placed = {p: classify(kb, p, s, leads, {}) for p, s in cands.items()}
    assert placed["9WF2"].placement != "related"
    assert {p: (x.placement, x.owner_pdb) for p, x in placed.items() if p != "9WF2"} == \
        {p: ("related", "9WF2") for p in ("9WF3", "9WFL", "9WFO")}


def test_note_and_evidence_say_why_an_old_structure_is_there():
    sibs = {"9WF2": {"paper": KEY, "with": ["9WFO"], "released": "2026-04-22"},
            "9WF3": {"paper": KEY, "with": ["9WFO", "9WFP"], "released": "2026-04-01"}}
    batch = lambda pdb, owner: {"pdb": pdb, "ownerPdb": owner, "rule": "same-paper-batch",
                                "evidence": [f"Same primary citation as {owner}, released in the same batch", "x"]}
    props = [{"pdb": "9WF2", "ownerPdb": "9BUM", "rule": "same-protein", "evidence": ["Same protein(s)"]},
             batch("9WF3", "9WF2"), batch("9WFO", "9WF2"), batch("9ABD", "9ABC")]
    mark_siblings(props, sibs, WINDOW)
    assert props[0]["sibling"] == {"of": DOI, "with": ["9WFO"], "when": "before",
                                   "reason": "released before the 26-week window; same paper as 9WFO"}
    assert props[0]["evidence"] == ["Same protein(s)"]
    assert props[1]["sibling"]["reason"].endswith("same paper as 9WFO and 1 more")
    assert props[1]["evidence"] == ["Same primary citation as 9WF2; both were released before the window", "x"]
    assert "sibling" not in props[2]
    assert props[2]["evidence"][0] == "Same primary citation as 9WF2 (released 2026-04-22, before the window)"
    assert props[3]["evidence"][0].endswith("released in the same batch")
    assert siblings.note({"paper": "pmid:123", "with": ["9BBB"]}, ["2026-01-01", "2026-01-29"]) == \
        {"of": "pmid:123", "with": ["9BBB"], "when": "before",
         "reason": "released before the 4-week window; same paper as 9BBB"}


# ---------------------------------------------------------------- the weekly run's step

def test_window_candidates_use_the_fresh_citation_of_tracked_structures(monkeypatch):
    monkeypatch.setattr(rcsb, "entry_summaries", lambda ids: {
        "9WFO": {"rcsb_primary_citation": PAPER}, "9TBP": {"rcsb_primary_citation": {"journal_abbrev": "To be published"}},
        "9OUT": {"rcsb_primary_citation": {"journal_abbrev": "To be published"}}, "9NIL": None})
    got = automation.window_candidates(["9NIL", "9OUT", "9TBP", "9WFO"], {"9OUT": dict(PAPER, pdbx_database_id_DOI="10.1/out")})
    assert [(c["pdb"], c["doi"], c["awaiting"]) for c in got] == \
        [("9NIL", None, True), ("9OUT", "10.1/out", False), ("9TBP", None, True), ("9WFO", DOI, False)]


def test_find_siblings_scores_them_and_applies_the_keeping_rule():
    kb = KB(entries=[], by_pdb={"9INM": []})
    seen = {}

    def discover(cands, window, have, listed, decided):
        seen.update(have=have, listed=listed, decided=decided, window=window)
        return ({p: {"paper": KEY, "with": ["9WFO"], "released": "2026-04-01"} for p in ("9WF2", "9WF3", "9SOL")},
                {"papers": 1, "found": 3, "papersWithSiblings": 1, "seconds": 0.0})

    score = lambda ids: {"9WF2": {"p": 0.7, "why": ["a"]}, "9WF3": {"p": 0.05, "why": ["b"]}, "9SOL": {"p": 0.02, "why": []}}
    cores = lambda ids: {"9WFO": {"P38435"}, "9WF2": {"P38435"}, "9WF3": {"P38435"}, "9SOL": {"P00740"}}
    sibs, mem, counts = automation.find_siblings([{"pdb": "9WFO", "doi": DOI}], WINDOW, {"9WFO"}, kb, {"9MEN": "k"}, {"9DEC"},
                                                 discover=discover, score=score, cores=cores)
    assert seen == {"have": {"9WFO"}, "listed": {"9INM", "9MEN"}, "decided": {"9DEC"}, "window": WINDOW}
    assert {p: s["kept"] for p, s in sibs.items()} == {"9WF2": "score", "9WF3": "same protein as 9WFO"}
    assert set(mem) == {"9WF2", "9WF3"}
    assert counts["kept"] == 2 and counts["belowFloor"] == 1 and counts["keptBySameProtein"] == 1


# ---------------------------------------------------------------- the in-place command

def _cand(pdb, membrane, **more):
    return {"pdb": pdb, "title": f"t {pdb}", "released": "2026-07-15", "doi": None, "pmid": None, "uniprot": ["P1"],
            "placement": "parent", "subgroup": "G03-S001", "subgroupName": "X", "group": "G", "ownerPdb": None,
            "ownerKey": None, "confidence": 0.5, "rule": "new-protein", "evidence": ["e"], "alternatives": [],
            "membrane": membrane, "membraneWhy": ["w"], "mp": {"name": f"n {pdb}", "description": "cryo-EM structure"},
            **more}


def _fake_run(monkeypatch, sibs, mem):
    """automation.add_siblings with the PDB, mpstruc and the listing drafts faked."""
    monkeypatch.setattr(automation, "_full_summaries", lambda db: {})
    monkeypatch.setattr(automation, "build_kb", lambda db, S: KB(entries=[]))
    monkeypatch.setattr(automation, "mentioned_in_descriptions", lambda db: {})
    monkeypatch.setattr(automation, "strand_rules", lambda kb: {})
    monkeypatch.setattr(automation, "find_siblings",
                        lambda cands, window, have, kb, mentioned, decided, db, discover=None, fallback=None:
                        (sibs, mem, {"papers": 1, "found": len(sibs), "seconds": 0.1, "belowFloor": 0}))

    def propose(kb, ids, conf, model, keep=None):
        lead = min(ids)
        return [{"pdb": p, "title": f"t {p}", "released": "2026-04-01", "doi": DOI, "pmid": "41888148", "uniprot": ["P38435"],
                 "sig": {}, "taxa": [], "likelihood": 0.5, "reasons": [],
                 "placement": "daughter" if p == lead else "related", "subgroup": "G03-S157", "subgroupName": "Carboxylases",
                 "group": "G", "ownerPdb": "9BUM" if p == lead else lead, "ownerKey": "K" if p == lead else None,
                 "confidence": 0.7, "rule": "same-protein" if p == lead else "same-paper-batch",
                 "evidence": ["Same protein(s)"] if p == lead else [f"Same primary citation as {lead}, released in the same batch"],
                 "alternatives": []} for p in ids], {}

    def add_listings(props, db=None, summaries=None):
        kids = {}
        for p in props:
            if p["placement"] == "related":
                kids.setdefault(p["ownerPdb"], []).append(p["pdb"])
        for p in props:
            p["mp"] = {"name": "GGCX", "description": "cryo-EM structure" + "".join(f"<br /> x: {k}" for k in sorted(kids.get(p["pdb"], [])))}
            if p["placement"] == "related":
                p["mp"]["relatedLine"] = f"x: {p['pdb']}"

    monkeypatch.setattr(automation, "propose", propose)
    monkeypatch.setattr(automation, "add_listings", add_listings)


def test_add_siblings_groups_the_study_and_leaves_every_other_candidate_alone(tmp_path, monkeypatch):
    wfo = _cand("9WFO", 0.7, doi=DOI, pmid="41888148", uniprot=["P38435"], placement="daughter", ownerPdb="9BUM",
                ownerKey="K", rule="same-protein", publishedOn="2026-11-12", tmEvidence=[{"uniprot": "P38435"}],
                mp={"name": "GGCX with bound proP1", "description": "cryo-EM structure", "ligandNotes": [{"id": "CLR"}]})
    doc = {"generated": "2026-11-13T10:00+00:00", "window": WINDOW, "stats": {"released": 3, "candidates": 3},
           "candidates": [_cand("9AAA", 0.9, odd={"é": 1}), wfo, _cand("9ZZZ", 0.2, awaiting={"since": "2026-11-01"})],
           "inMpstruc": [{"pdb": "9INM"}]}
    path = tmp_path / "candidates.json"
    path.write_text(json.dumps(doc, separators=(",", ":")), encoding="utf-8")
    (tmp_path / "curation.json").write_text(json.dumps({"entries": {"9DEC": {"status": "approved"}}}))
    sibs = {p: {"paper": KEY, "with": ["9WFO"], "released": d, "kept": "score"}
            for p, d in (("9WF2", "2026-04-22"), ("9WF3", "2026-04-01"))}
    mem = {"9WF2": {"p": 0.7, "why": ["9 segments"], "tmTable": [{"uniprot": "P38435"}]},
           "9WF3": {"p": 0.8, "why": ["9 segments"], "tmTable": [{"uniprot": "P38435"}], "barrel": {"count": 12}}}
    _fake_run(monkeypatch, sibs, mem)
    res = automation.add_siblings(path, db=object())
    text = path.read_text(encoding="utf-8")
    new = json.loads(text)
    assert text == json.dumps(new, separators=(",", ":"))  # still the weekly run's compact format
    by = {c["pdb"]: c for c in new["candidates"]}
    assert [c["pdb"] for c in new["candidates"]] == ["9AAA", "9WF3", "9WF2", "9WFO", "9ZZZ"]  # most likely first, then ID
    assert by["9AAA"] == doc["candidates"][0] and by["9ZZZ"] == doc["candidates"][2]
    assert new["inMpstruc"] == doc["inMpstruc"] and new["generated"] == doc["generated"]
    # The lead is the lowest ID, an older sibling; its description lists every related line.
    lead = by["9WF2"]
    assert (lead["placement"], lead["ownerPdb"]) == ("daughter", "9BUM")
    assert lead["mp"]["description"] == "cryo-EM structure<br /> x: 9WF3<br /> x: 9WFO"
    assert lead["sibling"] == {"of": DOI, "with": ["9WFO"], "when": "before",
                               "reason": "released before the 26-week window; same paper as 9WFO"}
    assert lead["membrane"] == 0.7 and lead["membraneWhy"] == ["9 segments"] and lead["tmEvidence"] == [{"uniprot": "P38435"}]
    assert not {"sig", "taxa", "likelihood", "reasons"} & set(lead)
    assert by["9WF3"]["barrelStrands"] == {"count": 12}
    assert by["9WF3"]["evidence"] == ["Same primary citation as 9WF2; both were released before the window"]
    # The candidate that was there moves under the new lead and keeps what the placement does not own.
    old = by["9WFO"]
    assert (old["placement"], old["ownerPdb"], old["rule"]) == ("related", "9WF2", "same-paper-batch")
    assert old["evidence"] == ["Same primary citation as 9WF2 (released 2026-04-22, before the window)"]
    assert old["publishedOn"] == "2026-11-12" and old["released"] == "2026-07-15" and old["title"] == "t 9WFO"
    assert old["mp"]["relatedLine"] == "x: 9WFO" and old["mp"]["ligandNotes"] == [{"id": "CLR"}]
    assert "sibling" not in old
    assert new["stats"] == {"released": 3, "candidates": 5, siblings.STAT: 2}
    assert [a["pdb"] for a in res["added"]] == ["9WF2", "9WF3"]
    assert res["studies"] == [{"paper": KEY, "before": ["9WFO daughter of 9BUM"],
                               "after": ["9WF2 daughter of 9BUM", "9WF3 related of 9WF2", "9WFO related of 9WF2"]}]


def test_add_siblings_writes_nothing_when_there_is_nothing_to_add(tmp_path, monkeypatch):
    doc = {"window": WINDOW, "stats": {"candidates": 1}, "candidates": [_cand("9AAA", 0.9)]}
    path = tmp_path / "candidates.json"
    path.write_text(json.dumps(doc, indent=1))  # not the compact format: left exactly as it is
    before = path.read_text()
    _fake_run(monkeypatch, {}, {})
    res = automation.add_siblings(path, db=object())
    assert res["added"] == [] and path.read_text() == before


# ---------------------------------------------------------------- review fixes and the user's rule

def test_lead_is_the_first_deposited_not_the_string_minimum_for_reissued_ids():
    from mpstruc.automation import lead_key
    s = lambda rel, dep: {"released": rel, "deposited": dep}
    # 24UG (deposited 2026) sorts below 7FBQ (2021) as a string but is the newer entry.
    assert min(["24UG", "7FBQ"], key=lambda p: lead_key(p, {"24UG": s("2026-06-17", "2026-03-20"),
                                                           "7FBQ": s("2022-07-13", "2021-07-12")}[p])) == "7FBQ"
    d = {"28PN": s("2026-05-06", "2026-02-12"), "28PP": s("2026-05-06", "2026-02-12"), "24UG": s("2026-06-17", "2026-03-20"),
         "9ST9": s("2026-05-06", "2025-09-26"), "9WYK": s("2026-05-20", "2025-09-27")}
    # Classic IDs first in string order, then the reissued ones by deposit date.
    assert sorted(d, key=lambda p: lead_key(p, d[p])) == ["9ST9", "9WYK", "28PN", "28PP", "24UG"]
    # An old classic 1xxx-3xxx ID is a classic ID.
    assert lead_key("1ABC", s("1999-01-01", "1998-01-01"))[0] == 0


def test_batch_leads_use_lead_order_and_same_name_strain_variants():
    kb = KB(entries=[])
    nfeob = lambda u, rel, dep: {**_s([u], [KEY], rel), "deposited": dep,
                                 "uniprot_names": {u: "ferrous iron transport protein b|vibrio cholerae"}}
    cands = {"9PS9": nfeob("A0A655NVH2", "2026-06-03", "2025-06-01"), "9PSD": nfeob("A0A655NVH2", "2026-06-03", "2025-06-01"),
             "9D8B": nfeob("A0A0X1KYD0", "2025-09-03", "2024-08-01"), "9D8D": nfeob("A0A0X1KYD0", "2025-09-03", "2024-08-01"),
             "28PN": {**_s(["P9"], [KEY], "2026-05-06"), "deposited": "2026-02-12"},
             "9ST9": {**_s(["P9"], [KEY], "2026-05-06"), "deposited": "2025-09-26"}}
    leads = batch_leads(cands)
    placed = {p: classify(kb, p, s, leads, {}) for p, s in cands.items()}
    assert {p: x.owner_pdb for p, x in placed.items() if x.placement == "related"} == \
        {"9D8D": "9D8B", "9PS9": "9D8B", "9PSD": "9D8B", "28PN": "9ST9"}


def test_protein_ids_ignore_generic_and_fusion_names():
    from mpstruc.automation import protein_ids
    s = {"core": {"U1", "U2", "U3"}, "uniprot_names": {"U1": "porin|pseudomonas aeruginosa",
                                                       "U2": "soluble cytochrome b562,receptor x|escherichia coli",
                                                       "U3": "ferrous iron transport protein b|"}}
    assert protein_ids(s) == {"U1", "U2", "U3"}


def test_select_takes_entries_released_after_the_window_too():
    found = {KEY: ["9WF2", "9NEW", "9LAT", "9LOW"]}
    released = {"9WF2": "2026-04-22", "9NEW": "2026-06-01", "9LAT": "2026-11-18", "9LOW": "2026-07-01"}
    sibs = siblings.select({KEY: ["9WFO"]}, found, released, WINDOW, have={"9WFO", "9NEW"}, listed=set(), decided=set())
    # 9NEW is a window release the run has; 9LOW one it judged (under the floor): not taken again.
    assert sorted(sibs) == ["9LAT", "9WF2"]
    lat = siblings.note({"paper": KEY, "with": ["9WFO"], "released": "2026-11-18"}, WINDOW)
    assert lat["when"] == "after" and lat["reason"].startswith("released after the 26-week window")
    low = siblings.note({"paper": KEY, "with": ["9WFO"], "released": "2026-03-01", "kept": "same protein as 9WFO"}, WINDOW)
    assert low["when"] == "before" and low["sameProteinAs"] == "9WFO"
    assert low["reason"].endswith("under 10% likely a membrane protein, kept because it is the same protein as 9WFO")


def test_failed_lookup_falls_back_to_the_expired_cache_then_to_the_previous_run(tmp_path, capsys):
    cache = tmp_path / "paper_entries.json"
    other = "doi:10.1/other"
    cache.write_text(json.dumps({KEY: {"at": 0.0, "ids": ["9WF2", "9WFO"]}}))
    tries = []

    def down(q):
        tries.append(q)
        raise TimeoutError("timed out")

    now = 30 * 86400.0  # the cache entry is long expired
    found, got = siblings.find_entries([KEY, other], down, None, cache, now=now, pause=0,
                                       fallback={other: ["9aaa", "9BBB"]})
    assert len(tries) == 2  # asked twice before giving up
    assert found == {KEY: ["9WF2", "9WFO"], other: ["9AAA", "9BBB"]} and got == {}
    assert json.loads(cache.read_text())[KEY]["at"] == 0.0  # not renewed: asked again next run
    assert "asking again" in capsys.readouterr().out
    # A batch that fails once and then answers is used.
    calls = iter([TimeoutError("x"), ["9WF2"]])

    def flaky(q):
        r = next(calls)
        if isinstance(r, Exception):
            raise r
        return r
    found, _ = siblings.find_entries([KEY], flaky, lambda ids: {"9WF2": PAPER}, cache, now=now, pause=0)
    assert found == {KEY: ["9WF2"]}


def test_the_cache_life_is_well_under_the_weekly_cadence():
    assert siblings.MAX_AGE_S <= 4 * 86400


def test_known_entries_are_every_structure_of_each_paper_in_the_file():
    cands = [{"pdb": "9wfo", "doi": DOI}, {"pdb": "9WF2", "doi": DOI.upper(), "sibling": {"of": DOI}},
             {"pdb": "9BBB", "doi": None, "pmid": "123"}, {"pdb": "9CCC"}]
    assert siblings.known_entries(cands) == {KEY: ["9WF2", "9WFO"], "pmid:123": ["9BBB"]}


def test_papers_of_window_releases_already_in_mpstruc_bring_their_other_structures_in(tmp_path, monkeypatch):
    # 9WFO was added to mpstruc on its own; its paper's other structures (9WF2, ...) must still be found.
    doc = {"generated": "2026-11-19T10:00+00:00", "window": WINDOW, "stats": {"candidates": 1},
           "candidates": [_cand("9AAA", 0.9)], "inMpstruc": [{"pdb": "9WFO", "as": "daughter"}]}
    path = tmp_path / "candidates.json"
    path.write_text(json.dumps(doc, separators=(",", ":")), encoding="utf-8")
    seen = {}
    _fake_run(monkeypatch, {}, {})
    monkeypatch.setattr(automation, "window_candidates",
                        lambda ids, fresh=None: [{"pdb": p, "doi": DOI, "pmid": "41888148", "awaiting": False} for p in ids])

    def find(cands, window, have, kb, mentioned, decided, db, discover=None, fallback=None):
        seen["cands"], seen["have"] = [c["pdb"] for c in cands], set(have)
        return {}, {}, {"papers": 1, "found": 0, "seconds": 0.1, "belowFloor": 0}
    monkeypatch.setattr(automation, "find_siblings", find)
    automation.add_siblings(path, db=object())
    assert seen["cands"] == ["9AAA", "9WFO"] and "9WFO" in seen["have"]
