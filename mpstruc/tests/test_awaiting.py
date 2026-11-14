from mpstruc import awaiting
from mpstruc.awaiting import ACCEPTED, PREPRINT_ONLY, TO_BE_PUBLISHED, citation_status, update


def test_citation_status_unpublished_any_case():
    for j in ("To Be Published", "to be published", "TO BE PUBLISHED", "To be Published.", "Not published"):
        assert citation_status({"journal_abbrev": j}) == TO_BE_PUBLISHED
    assert citation_status(None) == TO_BE_PUBLISHED
    assert citation_status({"journal_abbrev": None, "pdbx_database_id_DOI": None}) == TO_BE_PUBLISHED


def test_citation_status_preprints():
    for j in ("Biorxiv", "bioRxiv", "BioRXiv", "medRxiv", "Res Sq", "Research Square", "Chemrxiv", "arXiv", "SSRN", "Preprints"):
        assert citation_status({"journal_abbrev": j}) == PREPRINT_ONLY, j
    assert citation_status({"journal_abbrev": None, "pdbx_database_id_DOI": "10.1101/2025.03.04.641234"}) == PREPRINT_ONLY
    # bioRxiv's newer DOI prefix.
    assert citation_status({"journal_abbrev": None, "pdbx_database_id_DOI": "10.64898/2026.06.13.732081"}) == PREPRINT_ONLY


def test_citation_status_published():
    assert citation_status({"journal_abbrev": "Nat Commun", "pdbx_database_id_DOI": "10.1038/x"}) is None
    assert citation_status({"journal_abbrev": "Acta Pharmacol.Sin.", "journal_volume": "47"}) is None
    # Cold Spring Harbor journals share bioRxiv's DOI prefix.
    assert citation_status({"journal_abbrev": "Genes Dev.", "pdbx_database_id_DOI": "10.1101/gad.1234"}) is None
    assert citation_status({"journal_abbrev": None, "pdbx_database_id_DOI": "10.1038/s41586-026-1"}) is None


TBP = {"journal_abbrev": "To Be Published"}
PAPER = {"journal_abbrev": "Nature", "pdbx_database_id_DOI": "10.1038/x"}


def _c(pdb, released="2026-11-04"):
    return {"pdb": pdb, "released": released}


def test_first_seen_awaiting_and_since_kept():
    cands = [_c("22XQ"), _c("8ABC")]
    rec, stats = update(cands, {"22XQ": TBP, "8ABC": PAPER}, awaiting.empty(), "2026-11-05")
    assert cands[0]["awaiting"] == {"since": "2026-11-05", "reason": TO_BE_PUBLISHED}
    assert "awaiting" not in cands[1] and "publishedOn" not in cands[1]
    assert rec["awaiting"] == {"22XQ": {"since": "2026-11-05", "released": "2026-11-04", "reason": TO_BE_PUBLISHED}}
    assert stats == {"awaiting": 1, "newly_published": 0}
    # A week later it is still waiting: 'since' stays the first date.
    cands = [_c("22XQ")]
    rec, _ = update(cands, {"22XQ": TBP}, rec, "2026-11-12")
    assert cands[0]["awaiting"]["since"] == "2026-11-05"


def test_newly_published_gets_date_that_persists_then_expires():
    rec = {"awaiting": {"22XQ": {"since": "2026-11-05", "released": "2026-11-04", "reason": TO_BE_PUBLISHED}}, "published": {}}
    cands = [_c("22XQ")]
    rec, stats = update(cands, {"22XQ": PAPER}, rec, "2026-11-12")
    assert "awaiting" not in cands[0] and cands[0]["publishedOn"] == "2026-11-12"
    assert stats["newly_published"] == 1 and rec["awaiting"] == {}
    assert rec["published"]["22XQ"]["on"] == "2026-11-12"
    cands = [_c("22XQ")]
    rec, stats = update(cands, {"22XQ": PAPER}, rec, "2026-11-19")
    assert cands[0]["publishedOn"] == "2026-11-12" and stats["newly_published"] == 0
    cands = [_c("22XQ")]
    rec, _ = update(cands, {"22XQ": PAPER}, rec, "2027-06-01")
    assert "publishedOn" not in cands[0] and rec["published"] == {}


def test_waiting_over_two_years_is_given_up():
    rec = {"awaiting": {"22XQ": {"since": "2026-11-05", "released": "2026-11-04", "reason": TO_BE_PUBLISHED},
                        "7OLD": {"since": "2025-01-01", "released": "2024-10-01", "reason": TO_BE_PUBLISHED},
                        "7NRD": {"since": "2024-06-01", "reason": PREPRINT_ONLY}},
           "published": {"8PUB": {"on": "2026-10-01", "released": "2024-01-01"}}}
    kept, gone = awaiting.drop_stale(rec, "2026-11-12")
    assert gone == ["7NRD", "7OLD"]
    assert set(kept["awaiting"]) == {"22XQ"} and kept["published"] == rec["published"]
    assert awaiting.with_tracked([], kept) == ["22XQ", "8PUB"]
    # Just under two years: still tracked.
    _, gone = awaiting.drop_stale(rec, "2026-09-30")
    assert gone == ["7NRD"]


def test_forget_changed_leaves_expired_caches_alone_and_keeps_age(tmp_path, monkeypatch):
    import json, os, time
    from mpstruc import rcsb
    monkeypatch.setattr(rcsb, "CACHE_DIR", tmp_path)
    det, ent = tmp_path / "entry_details_v2.json", tmp_path / "entry_entities.json"
    det.write_text(json.dumps({"22XQ": {"rcsb_primary_citation": TBP}, "9OLD": {}}))
    ent.write_text(json.dumps({"22XQ": {}, "9OLD": {}}))
    three_days = time.time() - 3 * 86400
    os.utime(det, (three_days, three_days))
    os.utime(ent, (three_days, three_days))
    assert awaiting.forget_changed({"22XQ": PAPER}) == ["22XQ"]
    assert "22XQ" not in json.loads(det.read_text()) and "22XQ" not in json.loads(ent.read_text())
    assert abs(det.stat().st_mtime - three_days) < 1 and abs(ent.stat().st_mtime - three_days) < 1
    # Past entry_details' 7-day life: it is not rewritten (or made to look new).
    ten_days = time.time() - 10 * 86400
    det.write_text(json.dumps({"22XQ": {"rcsb_primary_citation": TBP}}))
    os.utime(det, (ten_days, ten_days))
    assert awaiting.forget_changed({"22XQ": PAPER}) == []
    assert abs(det.stat().st_mtime - ten_days) < 1


def test_citations_for_forgets_cached_entries_whose_paper_came_out(monkeypatch):
    from mpstruc import rcsb
    monkeypatch.setattr(rcsb, "entry_details", lambda ids: {"9NEW": {"rcsb_primary_citation": TBP}})
    monkeypatch.setattr(awaiting, "fresh_citations", lambda ids: {"9NEW": PAPER})
    forgot = []
    monkeypatch.setattr(awaiting, "forget_changed", lambda fresh: forgot.append(sorted(fresh)) or sorted(fresh))
    assert awaiting.citations_for([{"pdb": "9NEW"}], {}) == {"9NEW": PAPER}
    assert forgot == [["9NEW"]]


def test_unknown_citation_keeps_previous_state_and_decided_are_dropped():
    rec = {"awaiting": {"22XQ": {"since": "2026-11-05", "reason": PREPRINT_ONLY}, "22XR": {"since": "2026-11-05"}},
           "published": {}}
    cands = [_c("22XQ")]  # 22XR was approved or rejected, so it is no longer a candidate
    rec, _ = update(cands, {}, rec, "2026-11-12")
    assert cands[0]["awaiting"] == {"since": "2026-11-05", "reason": PREPRINT_ONLY}
    assert set(rec["awaiting"]) == {"22XQ"}


def test_tracked_ids_rejoin_the_run(tmp_path):
    rec = {"updated": None, "awaiting": {"1OLD": {"since": "2025-01-01"}}, "published": {"2OLD": {"on": "2026-11-01"}}}
    assert awaiting.with_tracked(["9NEW", "1OLD"], rec) == ["9NEW", "1OLD", "2OLD"]
    p = tmp_path / "awaiting.json"
    awaiting.save(p, rec)
    assert awaiting.load(p)["awaiting"] == rec["awaiting"]
    assert awaiting.load(tmp_path / "missing.json") == awaiting.empty()
    p.write_text('{"22xq": {"since": "2026-11-05", "released": "2026-11-04"}}')
    assert awaiting.load(p)["awaiting"] == {"22XQ": {"since": "2026-11-05", "released": "2026-11-04"}}


def test_a_journal_name_alone_is_accepted_not_published():
    # 9I1Z: the PDB swapped 'To be published' for 'Nat Commun' 2026, with no DOI, PubMed ID or pages.
    cit = {"title": "Membrane Protein Transporter", "journal_abbrev": "Nat Commun", "year": 2026}
    assert citation_status(cit) == ACCEPTED
    # Found by its title (the PDB had not added the DOI yet): published.
    assert citation_status({**cit, "found": {"pmid": "42919912", "doi": "10.1038/x", "source": "pubmed"}}) is None


def test_resolve_titles_and_fill_found(monkeypatch):
    from mpstruc import pubmed
    papers = {"Sustained antidepressant-like effects of psychedelics": {"pmid": "42919912", "doi": "10.1038/s41593-026-02497-5",
                                                                       "source": "pubmed"}}
    monkeypatch.setattr(pubmed, "find_paper", lambda title, year: papers.get(title))
    cits = {"9WZO": {"title": "Sustained antidepressant-like effects of psychedelics", "journal_abbrev": "Nat.Neurosci.",
                     "year": 2026},
            "9I1Z": {"title": "Membrane Protein Transporter", "journal_abbrev": "Nat Commun", "year": 2026},
            "9AAA": {"title": "x", "journal_abbrev": "Nature", "pdbx_database_id_DOI": "10.1038/y"}}
    assert awaiting.resolve_titles(cits) == 1
    assert citation_status(cits["9WZO"]) is None and citation_status(cits["9I1Z"]) == ACCEPTED
    cands = [{"pdb": "9WZO", "pmid": None, "doi": None, "mp": {"bibliography": {"pubMedId": None, "doi": None}}}]
    awaiting.fill_found(cands, cits)
    assert cands[0]["pmid"] == "42919912" and cands[0]["mp"]["bibliography"]["doi"] == "10.1038/s41593-026-02497-5"


def test_find_paper_needs_the_exact_title_and_year(monkeypatch, tmp_path):
    from mpstruc import pubmed
    monkeypatch.setattr(pubmed, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(pubmed.time, "sleep", lambda s: None)
    monkeypatch.setattr(pubmed, "_fetch", lambda url: {"esearchresult": {"idlist": ["31637997"]}} if "esearch" in url
                        else {"message": {"items": []}})
    monkeypatch.setattr(pubmed, "summaries", lambda ids, max_age_s=0: {
        "31637997": {"title": "Membrane protein transporter structure.", "year": "2019", "doi": "10.1/old"}})
    # Same words, but a 2019 paper: a placeholder title never matches an old paper.
    assert pubmed.find_paper("Membrane Protein Transporter Structure", 2026) is None
    assert pubmed.find_paper("Membrane Protein", 2026) is None  # too short to look up
