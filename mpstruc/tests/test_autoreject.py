import copy
import json
import shutil
import subprocess

import pytest

from mpstruc import autoreject, awaiting
from mpstruc.autoreject import (AUTO_BY, AUTO_NAME, decided_ids, mark_candidates, track_unpublished,
                                update_curation, update_text)

AT = "2026-11-12T06:30:00.000Z"
TBP = {"since": "2026-11-05", "reason": awaiting.TO_BE_PUBLISHED}


def _cand(pdb, **kw):
    c = {"pdb": pdb, "title": f"structure {pdb}", "released": "2026-11-04", "organisms": ["Escherichia coli"],
         "res": 2.5, "method": "ELECTRON MICROSCOPY", "pmid": None, "doi": None,
         "mp": {"name": "permease", "bibliography": {"title": "x"}}, "placement": "related",
         "subgroup": "G03-S084", "ownerPdb": "22XP", "membrane": 0.95}
    c.update(kw)
    return c


def _rec(wait=(), pub=None):
    return {"updated": "2026-11-12", "awaiting": {p: {"since": "2026-11-05", "reason": awaiting.TO_BE_PUBLISHED}
                                                  for p in wait}, "published": pub or {}}


def _curator_reject(reason="Other", code=None):
    e = {"status": "rejected", "reason": reason, "by": "jaylinantwerp-beep", "byName": "jaylinantwerp-beep",
         "at": "2026-11-06T00:11:16.810Z", "title": "structure of a permease", "released": "2026-11-04"}
    if code:
        e["reasonCode"] = code
    return e


def _raw(entries=None, log=None):
    return {"version": 1, "entries": entries or {}, "log": log or []}


def test_waiting_suggestion_is_rejected_automatically():
    raw = _raw()
    cands = [_cand("22XQ", awaiting=TBP), _cand("8ABC"),
             _cand("9PRE", awaiting={"since": "2026-11-05", "reason": awaiting.PREPRINT_ONLY})]
    counts = update_curation(raw, cands, _rec(["22XQ", "9PRE"]), AT)
    assert counts == {"rejected": ["22XQ", "9PRE"], "unrejected": []}
    e = raw["entries"]["22XQ"]
    assert e["status"] == "rejected" and e["reasonCode"] == "unpublished" and e["auto"] is True
    assert e["reason"] == "Not published yet (PDB lists it as 'To be published')"
    assert raw["entries"]["9PRE"]["reason"] == "Not published yet (preprint only)"
    assert (e["by"], e["byName"], e["at"]) == (AUTO_BY, AUTO_NAME, AT)
    assert e["title"] == "structure 22XQ" and e["organism"] == "Escherichia coli" and e["suggestedOwner"] == "22XP"
    assert e["mp"] == cands[0]["mp"] and e["suggested"] == "related" and e["membrane"] == 0.95
    assert "8ABC" not in raw["entries"]
    line = raw["log"][0]
    assert line["action"] == "auto-reject" and line["pdb"] == "22XQ" and line["by"] == AUTO_BY
    assert line["detail"] == ("Rejected 22XQ automatically: not published yet (PDB lists it as 'To be published'). "
                              "It will come back to To decide when the paper is out.")


def test_published_rejections_come_back_auto_and_curator():
    ents = {"22XQ": {"status": "rejected", "reasonCode": "unpublished", "reason": "x", "by": AUTO_BY, "auto": True},
            "22XR": _curator_reject("Not published yet", "unpublished"),
            "22XS": _curator_reject("Not published yet", "unpublished")}
    raw = _raw(ents)
    cands = [_cand("22XQ", publishedOn="2026-11-12"), _cand("22XS", awaiting=TBP)]
    # 22XR is no longer a candidate (e.g. added to mpstruc) but the record saw its paper.
    rec = _rec(["22XS"], {"22XQ": {"on": "2026-11-12"}, "22XR": {"on": "2026-11-10"}})
    counts = update_curation(raw, cands, rec, AT)
    assert counts == {"rejected": [], "unrejected": ["22XQ", "22XR"]}
    assert list(raw["entries"]) == ["22XS"]
    assert [l["detail"] for l in raw["log"]] == ["Paper for 22XQ is now published (2026-11-12): back in To decide",
                                                "Paper for 22XR is now published (2026-11-10): back in To decide"]
    assert all(l["action"] == "auto-unreject" for l in raw["log"])


def test_never_touches_tombstones_approved_or_other_reasons():
    ents = {"1DEL": {"status": "deleted", "reasonCode": "unpublished", "reason": "x", "deletedFrom": "rejected"},
            "2APP": {"status": "approved", "placement": "parent"},
            "3OTH": _curator_reject("Other"),                       # old entry: no reasonCode
            "4DUP": _curator_reject("Duplicate", "duplicate"),
            "22XQ": _curator_reject("Other")}                       # migration: curator's, left as is
    raw = _raw(copy.deepcopy(ents))
    cands = [_cand(p, awaiting=TBP) for p in ("1DEL", "2APP", "3OTH", "4DUP", "22XQ")]
    pub = {p: {"on": "2026-11-12"} for p in ents}
    assert update_curation(raw, cands, _rec(list(ents)), AT) == {"rejected": [], "unrejected": []}
    cands = [_cand(p, publishedOn="2026-11-12") for p in ents]
    assert update_curation(raw, cands, _rec([], pub), AT) == {"rejected": [], "unrejected": []}
    assert raw == _raw(ents)


def test_idempotent_across_runs():
    raw = _raw({"2APP": {"status": "approved"}}, [{"at": "x", "by": "a", "action": "approve", "pdb": "2APP"}])
    cands = [_cand("22XQ", awaiting=TBP)]
    update_curation(raw, cands, _rec(["22XQ"]), AT)
    once = copy.deepcopy(raw)
    assert update_curation(raw, cands, _rec(["22XQ"]), "2026-11-19T06:30:00.000Z") == {"rejected": [], "unrejected": []}
    assert raw == once
    # Paper out: back once, and not rejected again next week.
    cands = [_cand("22XQ", publishedOn="2026-11-26")]
    rec = _rec([], {"22XQ": {"on": "2026-11-26"}})
    assert update_curation(raw, cands, rec, AT)["unrejected"] == ["22XQ"]
    after = copy.deepcopy(raw)
    assert update_curation(raw, cands, rec, AT) == {"rejected": [], "unrejected": []}
    assert raw == after


def test_still_waiting_or_unknown_stays_rejected():
    raw = _raw({"22XQ": _curator_reject("Not published yet", "unpublished")})
    assert update_curation(raw, [_cand("22XQ", awaiting=TBP)], _rec(["22XQ"]), AT)["unrejected"] == []
    # Not a candidate and nothing known about its paper.
    assert update_curation(raw, [], _rec(), AT)["unrejected"] == []
    assert "22XQ" in raw["entries"]


def test_curator_unarchive_is_respected():
    log = [{"by": AUTO_BY, "action": "auto-reject", "pdb": "22XQ"},
           {"by": "andyyjiang", "action": "unarchive", "pdb": "22XQ"}]
    raw = _raw({}, copy.deepcopy(log))
    assert update_curation(raw, [_cand("22XQ", awaiting=TBP)], _rec(["22XQ"]), AT)["rejected"] == []
    # An automatic rejection after the unarchive (e.g. after a new upload) is not blocked.
    raw = _raw({}, log + [{"by": AUTO_BY, "action": "auto-reject", "pdb": "22XQ"}])
    assert update_curation(raw, [_cand("22XQ", awaiting=TBP)], _rec(["22XQ"]), AT)["rejected"] == ["22XQ"]


def test_update_text_preserves_everything_else():
    doc = {"version": 1, "zzz": {"keep": ["me", 1.5]}, "entries": {
        "22XS": {"status": "approved", "by": "andyyjiang", "mpEdited": {"name": "Ünïcode"}, "edited": True},
        "23JK": _curator_reject("Other")}, "log": [{"at": "1", "action": "approve", "pdb": "22XS"}], "extra": None}
    text = json.dumps(doc, indent=1, ensure_ascii=False) + "\n"
    same, _, counts = update_text(text, [_cand("22XS", awaiting=TBP)], _rec(["22XS"]), AT)
    assert same is text and counts == {"rejected": [], "unrejected": []}
    new, raw, counts = update_text(text, [_cand("35YU", awaiting=TBP)], _rec(["35YU"]), AT)
    assert counts["rejected"] == ["35YU"]
    got = json.loads(new)
    assert list(got) == list(doc)
    assert list(got["entries"]) == ["22XS", "23JK", "35YU"]
    assert got["entries"]["22XS"] == doc["entries"]["22XS"] and got["entries"]["23JK"] == doc["entries"]["23JK"]
    assert got["zzz"] == doc["zzz"] and got["extra"] is None and got["log"][0] == doc["log"][0]
    # Unchanged lines stay byte-for-byte (the site's own format).
    prefix = text[:text.index('"released": "2026-11-04"')]
    assert new.startswith(prefix) and new.endswith("\n")
    assert "Ünïcode" in new
    with pytest.raises(ValueError):
        update_text('{"entries": []}', [], _rec(), AT)
    new, _, _ = update_text(None, [_cand("35YU", awaiting=TBP)], _rec(["35YU"]), AT)
    assert json.loads(new)["version"] == 1


def test_decided_and_tracking_keep_unpublished_rejections_checked():
    cur = _raw({"22XQ": _curator_reject("Not published yet", "unpublished"), "22XR": _curator_reject("Other"),
                "1DEL": {"status": "deleted", "reasonCode": "unpublished"}, "2APP": {"status": "approved"},
                "7OLD": dict(_curator_reject("Not published yet", "unpublished"), released="2023-01-01"),
                "9PRE": _curator_reject("Not published yet (preprint only)", "unpublished")})
    assert decided_ids(cur) == {"22XR", "1DEL", "2APP"}
    rec, added = track_unpublished(awaiting.empty(), cur, "2026-11-12")
    assert added == ["22XQ", "9PRE"]  # 7OLD released over two years ago: given up
    assert rec["awaiting"]["22XQ"] == {"since": "2026-11-06", "released": "2026-11-04",
                                       "reason": awaiting.TO_BE_PUBLISHED, "unconfirmed": True}
    assert rec["awaiting"]["9PRE"]["reason"] == awaiting.PREPRINT_ONLY
    # Already tracked: left as is.
    rec2, added = track_unpublished(rec, cur, "2026-11-12")
    assert added == [] and rec2 == rec


def test_awaiting_update_holds_rejected_structures():
    rec = {"awaiting": {"22XQ": {"since": "2026-11-05", "released": "2026-11-04", "reason": awaiting.TO_BE_PUBLISHED},
                        "22XR": {"since": "2026-11-05", "released": "2026-11-04", "reason": awaiting.TO_BE_PUBLISHED},
                        "22XT": {"since": "2026-11-05", "released": "2026-11-04", "reason": awaiting.TO_BE_PUBLISHED}},
           "published": {}}
    paper = {"journal_abbrev": "Nature", "pdbx_database_id_DOI": "10.1038/x"}
    # None of them is a candidate this week; 22XQ and 22XR are rejected as unpublished, 22XT is not.
    new, stats = awaiting.update([], {"22XR": paper}, rec, "2026-11-12", hold={"22XQ", "22XR"})
    assert set(new["awaiting"]) == {"22XQ"} and new["awaiting"]["22XQ"]["since"] == "2026-11-05"
    assert new["published"] == {"22XR": {"on": "2026-11-12", "released": "2026-11-04"}}
    assert stats == {"awaiting": 1, "newly_published": 1}
    # The rejection of 22XR is then lifted.
    raw = _raw({"22XR": _curator_reject("Not published yet", "unpublished")})
    assert update_curation(raw, [], new, AT)["unrejected"] == ["22XR"]


PAPER = {"journal_abbrev": "Nature", "pdbx_database_id_DOI": "10.1038/x", "page_first": None,
         "pdbx_database_id_PubMed": None}


def _week(raw, cands, cits, rec, today):
    """The order automation.main() runs things in, on pure data."""
    hold = set(autoreject.unpublished_rejections(raw))
    rec, _ = track_unpublished(rec, raw, today)
    rec, _ = awaiting.drop_stale(rec, today)
    rec, stats = awaiting.update(cands, cits, rec, today, hold=hold)
    counts = update_curation(raw, cands, rec, f"{today}T06:30:00.000Z")
    mark_candidates(cands, raw)
    return rec, stats, counts


def test_curator_rejection_against_a_journal_citation_sticks():
    # The PDB already names a journal (no page numbers yet); a curator rejects
    # it as "Not published yet". The run must not call it newly published.
    ent = dict(_curator_reject("Not published yet: in press, no page numbers", "unpublished"),
               at="2026-11-09T10:00:00.000Z")
    raw = _raw({"3CCC": ent}, [{"at": ent["at"], "by": ent["by"], "action": "reject", "pdb": "3CCC"}])
    rec = awaiting.empty()
    for day in ("2026-11-12", "2026-11-19"):
        cands = [_cand("3CCC")]
        rec, stats, counts = _week(raw, cands, {"3CCC": PAPER}, rec, day)
        assert counts == {"rejected": [], "unrejected": []}
        assert stats["newly_published"] == 0 and "publishedOn" not in cands[0]
        assert "3CCC" not in rec["published"] and raw["entries"]["3CCC"] == ent
    # Also when it is not among the week's suggestions.
    rec, stats, counts = _week(raw, [], {"3CCC": PAPER}, rec, "2026-11-26")
    assert counts["unrejected"] == [] and stats["newly_published"] == 0 and "3CCC" in raw["entries"]


def test_curator_rejection_comes_back_once_the_run_saw_it_waiting():
    ent = dict(_curator_reject("Not published yet", "unpublished"), at="2026-11-09T10:00:00.000Z")
    raw = _raw({"3CCC": ent})
    rec, _, counts = _week(raw, [_cand("3CCC")], {"3CCC": None}, awaiting.empty(), "2026-11-12")
    assert counts["unrejected"] == [] and "unconfirmed" not in rec["awaiting"]["3CCC"]
    # Unknown PDB data keeps the mark as it was (here: confirmed).
    cands = [_cand("3CCC")]
    rec, stats, counts = _week(raw, cands, {"3CCC": PAPER}, rec, "2026-11-19")
    assert counts["unrejected"] == ["3CCC"] and stats["newly_published"] == 1
    assert cands[0]["publishedOn"] == "2026-11-19"


def test_unconfirmed_mark_kept_while_pdb_data_unavailable():
    ent = dict(_curator_reject("Not published yet", "unpublished"), at="2026-11-09T10:00:00.000Z")
    raw = _raw({"3CCC": ent})
    rec, _, _ = _week(raw, [_cand("3CCC")], {}, awaiting.empty(), "2026-11-12")
    assert rec["awaiting"]["3CCC"]["unconfirmed"] is True
    rec, stats, counts = _week(raw, [_cand("3CCC")], {"3CCC": PAPER}, rec, "2026-11-19")
    assert counts["unrejected"] == [] and stats["newly_published"] == 0


def test_curator_rejecting_again_after_it_came_back_is_not_overruled():
    raw = _raw()
    rec = awaiting.empty()
    rec, _, counts = _week(raw, [_cand("9XYZ")], {"9XYZ": None}, rec, "2026-11-05")
    assert counts["rejected"] == ["9XYZ"]
    rec, _, counts = _week(raw, [_cand("9XYZ")], {"9XYZ": PAPER}, rec, "2026-11-12")
    assert counts["unrejected"] == ["9XYZ"] and "9XYZ" not in raw["entries"]
    # The curator wants page numbers first, and rejects it again.
    ent = dict(_curator_reject("Not published yet: no page numbers", "unpublished"), at="2026-11-13T09:00:00.000Z")
    raw["entries"]["9XYZ"] = ent
    raw["log"].append({"at": ent["at"], "by": ent["by"], "action": "reject", "pdb": "9XYZ"})
    for day in ("2026-11-19", "2026-11-26"):
        cands = [_cand("9XYZ")]
        rec, _, counts = _week(raw, cands, {"9XYZ": PAPER}, rec, day)
        assert counts == {"rejected": [], "unrejected": []} and raw["entries"]["9XYZ"] == ent
    # Same day as the automatic un-rejection: also left alone.
    raw2 = _raw({"9XYZ": dict(ent, at="2026-11-12T09:00:00.000Z")})
    assert update_curation(raw2, [_cand("9XYZ", publishedOn="2026-11-12")], rec, AT)["unrejected"] == []


@pytest.mark.parametrize("action", ["restore", "unarchive", "undo"])
def test_restore_from_deleted_or_undo_is_respected(action):
    raw = _raw()
    rec, _, counts = _week(raw, [_cand("4DDD")], {"4DDD": None}, awaiting.empty(), "2026-11-05")
    assert counts["rejected"] == ["4DDD"]
    if action == "restore":
        raw["entries"]["4DDD"] = {"status": "deleted", "reasonCode": "unpublished", "reason": "x",
                                  "by": "alice", "deletedFrom": "rejected"}
        raw["log"].append({"by": "alice", "action": "delete", "pdb": "4DDD"})
        rec, _, counts = _week(raw, [_cand("4DDD")], {"4DDD": None}, rec, "2026-11-12")
        assert counts == {"rejected": [], "unrejected": []}
    del raw["entries"]["4DDD"]
    raw["log"].append({"by": "alice", "action": action, "pdb": "4DDD"})
    for day in ("2026-11-19", "2026-11-26"):
        rec, _, counts = _week(raw, [_cand("4DDD")], {"4DDD": None}, rec, day)
        assert counts == {"rejected": [], "unrejected": []} and "4DDD" not in raw["entries"]


def test_undo_without_automatic_rejection_does_not_block():
    # A curator approved a waiting structure and undid it: the automation may still reject it.
    log = [{"by": "alice", "action": "approve", "pdb": "4DDD"}, {"by": "alice", "action": "undo", "pdb": "4DDD"}]
    raw = _raw({}, log)
    assert update_curation(raw, [_cand("4DDD", awaiting=TBP)], _rec(["4DDD"]), AT)["rejected"] == ["4DDD"]


def test_lone_surrogate_is_escaped_like_json_stringify(tmp_path):
    text = '{"version": 1, "entries": {"22XQ": {"status": "rejected", "reason": "Other: bad \\ud83d \\ud83d\\ude00"}}, "log": []}'
    new, raw, counts = update_text(text, [_cand("10ED", awaiting=TBP)], _rec(["10ED"]), AT)
    assert counts["rejected"] == ["10ED"]
    p = tmp_path / "curation.json"
    p.write_text(new, encoding="utf-8")
    assert '"Other: bad \\ud83d \U0001F600"' in p.read_text(encoding="utf-8")
    assert json.loads(p.read_text(encoding="utf-8"))["entries"]["22XQ"]["reason"] == raw["entries"]["22XQ"]["reason"]


def test_apply_updates_candidates_stats(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    gen = "2026-11-12T06:00:00Z"
    stats = {"rejected automatically: not published yet": 2, "back after publication": 0}
    (data / "candidates.json").write_text(json.dumps({"generated": gen, "stats": stats,
                                                      "candidates": [_cand("10ED", awaiting=TBP)]}))
    awaiting.save(data / "awaiting.json", _rec(["10ED"]))
    (data / "import_log.json").write_text(json.dumps({"runs": [{"at": gen, "autoRejected": 2}]}))
    (data / "curation.json").write_text(autoreject.dumps(_raw()))
    counts = autoreject.apply(data)
    assert counts["rejected"] == ["10ED"]
    out = json.loads((data / "candidates.json").read_text())
    assert out["stats"] == {"rejected automatically: not published yet": 1, "back after publication": 0}
    assert out["candidates"][0]["autoRejected"] is True
    assert json.loads((data / "import_log.json").read_text())["runs"][-1]["autoRejected"] == 1


def test_mark_candidates():
    raw = _raw({"22XQ": {"status": "rejected", "reasonCode": "unpublished", "auto": True},
                "22XR": _curator_reject("Not published yet", "unpublished")})
    cands = [_cand("22XQ"), _cand("22XR"), _cand("8ABC", autoRejected=True)]
    mark_candidates(cands, raw)
    assert [c.get("autoRejected") for c in cands] == [True, None, None]


def test_commit_message():
    assert autoreject.commit_message("2026-11-12", 2, 1) == \
        "Weekly mpstruc update (2026-11-12): 2 auto-rejected as not published yet, 1 back after publication"


# ------------------------------------------------------------ push and retry

def _git(cwd, *args):
    r = subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *args], cwd=cwd,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


@pytest.fixture
def repos(tmp_path):
    if not shutil.which("git"):
        pytest.skip("git not installed")
    origin = tmp_path / "origin.git"
    _git(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    seed = tmp_path / "seed"
    _git(tmp_path, "clone", "-q", str(origin), str(seed))
    (seed / "data").mkdir()
    cur = _raw({"23JK": _curator_reject("Other")}, [{"at": "1", "by": "j", "action": "reject", "pdb": "23JK"}])
    (seed / "data" / "curation.json").write_text(autoreject.dumps(cur))
    (seed / "data" / "candidates.json").write_text("{}")
    _git(seed, "add", ".")
    _git(seed, "commit", "-q", "-m", "seed")
    _git(seed, "push", "-q", "origin", "HEAD:main")
    job = tmp_path / "job"
    _git(tmp_path, "clone", "-q", str(origin), str(job))
    for k, v in (("user.name", "bot"), ("user.email", "bot@x")):
        _git(job, "config", k, v)
    return origin, seed, job


FILES = ["data/candidates.json", "data/curation.json"]


def _job_apply(job, cands, rec):
    path = job / "data" / "curation.json"
    new, _, _ = update_text(path.read_text(), cands, rec, AT)
    path.write_text(new)


def test_publish_retries_after_curation_conflict(repos):
    origin, seed, job = repos
    cands, rec = [_cand("22XQ", awaiting=TBP)], _rec(["22XQ"])
    # The job's run.
    (job / "data" / "candidates.json").write_text('{"candidates": 1}')
    _job_apply(job, cands, rec)
    # Meanwhile a curator saves a decision through the site (appends a log line too).
    path = seed / "data" / "curation.json"
    cur = json.loads(path.read_text())
    cur["entries"]["35YU"] = {"status": "approved", "by": "andyyjiang"}
    cur["log"].append({"at": "2", "by": "andyyjiang", "action": "approve", "pdb": "35YU"})
    path.write_text(autoreject.dumps(cur))
    _git(seed, "commit", "-q", "-am", "Curator andyyjiang: approve 35YU")
    _git(seed, "push", "-q", "origin", "HEAD:main")
    calls = []

    def reapply():
        calls.append(1)
        _job_apply(job, cands, rec)

    assert autoreject.publish(job, "main", FILES, lambda: "Weekly", reapply) == "pushed"
    assert calls == [1]  # the rebase conflicted on the log; re-applied once
    check = origin.parent / "check"
    _git(origin.parent, "clone", "-q", str(origin), str(check))
    got = json.loads((check / "data" / "curation.json").read_text())
    assert set(got["entries"]) == {"23JK", "35YU", "22XQ"}
    assert [l["action"] for l in got["log"]] == ["reject", "approve", "auto-reject"]
    assert (check / "data" / "candidates.json").read_text() == '{"candidates": 1}'
    assert _git(check, "log", "-1", "--format=%s").strip() == "Weekly"


def test_publish_without_conflict_and_without_changes(repos):
    origin, seed, job = repos
    assert autoreject.publish(job, "main", FILES, lambda: "x", lambda: None) == "no changes"
    (job / "data" / "candidates.json").write_text('{"candidates": 2}')
    # Someone else changes another file meanwhile: a plain rebase is enough.
    (seed / "README").write_text("hi")
    _git(seed, "add", "README")
    _git(seed, "commit", "-q", "-m", "docs")
    _git(seed, "push", "-q", "origin", "HEAD:main")
    calls = []
    assert autoreject.publish(job, "main", FILES, lambda: "Weekly", lambda: calls.append(1)) == "pushed"
    assert calls == []
    assert _git(job, "log", "--format=%s", "-2").split("\n")[:2] == ["Weekly", "docs"]
