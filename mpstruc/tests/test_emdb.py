"""EMDB map IDs: the RCSB query (rcsb.emdb_ids) and the in-place patch of the site files (web.patch_links)."""
import json

import pytest

from mpstruc import rcsb, web


def _fake_post(answers, calls):
    def post(url, payload):
        ids = payload["variables"]["ids"]
        calls.append(list(ids))
        entries = [None if pid not in answers else
                   {"rcsb_id": pid, "rcsb_entry_container_identifiers": {"emdb_ids": answers[pid]}} for pid in ids]
        return {"data": {"entries": [e for e in entries if e is not None]}}
    return post


def test_emdb_ids_batches_validates_and_caches(tmp_path, monkeypatch):
    monkeypatch.setattr(rcsb, "CACHE_DIR", tmp_path)
    answers = {"7ABC": ["EMD-12345", "EMD-12346", "EMD-12345"], "1PRH": None, "2BAD": ["EMD-1", "<b>", "EMD-0502"]}
    calls = []
    got = rcsb.emdb_ids(["7abc", "1PRH", "2BAD", "9ZZZ"], batch=2, post=_fake_post(answers, calls))
    assert got == {"7abc": ["EMD-12345", "EMD-12346"], "1PRH": [], "2BAD": ["EMD-0502"], "9ZZZ": []}
    assert calls == [["1PRH", "2BAD"], ["7ABC", "9ZZZ"]]  # 2 per query
    cached = json.loads((tmp_path / "emdb_ids.json").read_text())
    assert cached["ids"]["1PRH"] == [] and cached["ids"]["9ZZZ"] == []  # none is cached too
    # A second call asks only for what is new.
    calls.clear()
    rcsb.emdb_ids(["7ABC", "3NEW"], post=_fake_post({"3NEW": ["EMD-40000"]}, calls))
    assert calls == [["3NEW"]]
    # Once the cache is older than max_age_s everything is asked for again.
    calls.clear()
    rcsb.emdb_ids(["7ABC"], max_age_s=-1, post=_fake_post(answers, calls))
    assert calls == [["7ABC"]]


def test_pdb_info_has_emdb_only_when_there_are_maps():
    s = {"struct": {"title": "T"}, "exptl": [{"method": "ELECTRON MICROSCOPY"}],
         "rcsb_entry_info": {"resolution_combined": [3.1]}, "rcsb_accession_info": {"initial_release_date": "2020-01-01T00:00:00Z"}}
    assert web._pdb_info(s, ["EMD-1234"])["emdb"] == ["EMD-1234"]
    assert "emdb" not in web._pdb_info(s, [])
    assert "emdb" not in web._pdb_info(s)
    assert list(web._pdb_info(s, ["EMD-1234"]))[-1] == "emdb"


def _site(tmp_path):
    mp = {"meta": {"source": "Ö"}, "groups": [{"k": "G1", "subgroups": [{"k": "S1", "parents": [
        {"k": "P1", "pdb": "7ABC", "name": "α channel", "x": {"title": "cryo", "method": "Cryo-EM"},
         "d": [{"k": "D1", "pdb": "1prh", "x": {"title": "xray", "method": "X-ray Diffraction", "emdb": ["EMD-9999"]}},
               {"k": "D2", "pdb": "5OLD"}]}]}]}],
          "relPdb": {"6REL": {"title": "rel", "method": "Cryo-EM"}, "4XRY": {"title": "x"}}}
    cands = {"generated": "now", "candidates": [
        {"pdb": "8NEW", "title": "é", "mp": {"name": "n"}},
        {"pdb": "8SIB", "title": "s", "sibling": {"with": ["8NEW"]}},
        {"pdb": "8XRY", "title": "x"}]}
    mp_path, cand_path = tmp_path / "mpstruc.json", tmp_path / "candidates.json"
    # Written the way web.main and automation write them: compact, the first keeping non-ASCII as is.
    mp_path.write_text(json.dumps(mp, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    cand_path.write_text(json.dumps(cands, separators=(",", ":")), encoding="utf-8")
    return mp_path, cand_path


MAPS = {"7ABC": ["EMD-11111"], "1PRH": [], "6REL": ["EMD-22222", "EMD-22223"], "4XRY": [],
        "8NEW": ["EMD-33333"], "8SIB": ["EMD-44444", "javascript:x"], "8XRY": [], "5OLD": ["EMD-55555"]}


def test_patch_links_updates_every_record_in_place(tmp_path):
    mp_path, cand_path = _site(tmp_path)
    asked = []

    def fetch(ids):
        asked.append(ids)
        return {i: MAPS.get(i, []) for i in ids}

    counts = web.patch_links(mp_path, cand_path, fetch=fetch)
    assert asked == [sorted(["7ABC", "1PRH", "6REL", "4XRY", "8NEW", "8SIB", "8XRY"])]  # one batch call, D2 has no "x"
    mp_text = mp_path.read_text(encoding="utf-8")
    mp = json.loads(mp_text)
    p = mp["groups"][0]["subgroups"][0]["parents"][0]
    assert p["x"] == {"title": "cryo", "method": "Cryo-EM", "emdb": ["EMD-11111"]}
    assert "emdb" not in p["d"][0]["x"]  # stale map removed
    assert "x" not in p["d"][1]  # an entry the PDB does not know stays as it is
    assert mp["relPdb"]["6REL"]["emdb"] == ["EMD-22222", "EMD-22223"] and "emdb" not in mp["relPdb"]["4XRY"]
    assert mp["meta"] == {"source": "Ö"} and p["name"] == "α channel"
    assert mp_text == json.dumps(mp, ensure_ascii=False, separators=(",", ":"))  # same formatting
    cand_text = cand_path.read_text(encoding="utf-8")
    cands = json.loads(cand_text)
    assert [c.get("emdb") for c in cands["candidates"]] == [["EMD-33333"], ["EMD-44444"], None]
    assert cand_text == json.dumps(cands, separators=(",", ":")) and "\\u00e9" in cand_text
    assert counts["mpstruc.json"] == {"structures": 4, "withEmdb": 2, "recordsChanged": 3}
    assert counts["candidates.json"] == {"structures": 3, "withEmdb": 2, "recordsChanged": 2}


def test_patch_links_is_idempotent(tmp_path):
    mp_path, cand_path = _site(tmp_path)
    fetch = lambda ids: {i: MAPS.get(i, []) for i in ids}  # noqa: E731
    web.patch_links(mp_path, cand_path, fetch=fetch)
    before = (mp_path.stat().st_mtime_ns, mp_path.read_bytes(), cand_path.read_bytes())
    counts = web.patch_links(mp_path, cand_path, fetch=fetch)
    assert (mp_path.stat().st_mtime_ns, mp_path.read_bytes(), cand_path.read_bytes()) == before
    assert all(c["recordsChanged"] == 0 for c in counts.values())


def test_patch_links_keeps_indented_files_indented(tmp_path):
    mp_path, cand_path = _site(tmp_path)
    doc = json.loads(cand_path.read_text())
    cand_path.write_text(json.dumps(doc, indent=1) + "\n")
    web.patch_links(mp_path, cand_path, fetch=lambda ids: {i: MAPS.get(i, []) for i in ids})
    text = cand_path.read_text()
    assert text == json.dumps(json.loads(text), indent=1) + "\n"


def test_patch_links_refuses_unknown_formatting(tmp_path):
    mp_path, cand_path = _site(tmp_path)
    cand_path.write_text('{"candidates" :  []}')
    with pytest.raises(SystemExit):
        web.patch_links(mp_path, cand_path, fetch=lambda ids: {})
    assert cand_path.read_text() == '{"candidates" :  []}'


def test_links_command_runs_patch(tmp_path, monkeypatch, capsys):
    mp_path, cand_path = _site(tmp_path)
    monkeypatch.setattr(rcsb, "emdb_ids", lambda ids: {i: MAPS.get(i, []) for i in ids})
    web.main(["links", "--mpstruc", str(mp_path), "--candidates", str(cand_path)])
    out = capsys.readouterr().out
    assert "mpstruc.json: 2 of 4 structures have EMDB maps" in out
    assert json.loads(cand_path.read_text())["candidates"][0]["emdb"] == ["EMD-33333"]
