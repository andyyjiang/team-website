from types import SimpleNamespace

from mpstruc.membrane import _apply_isotonic, _isotonic, explain, features, predict


def _facts(**kw):
    f = {"pdb": "9NEW", "released": "2025-01-10", "uniprots": ["P1"], "tm_uniprots": ["P1"], "keywords": ["KW-0812"],
         "detergents": ["LMN"], "lipids": [], "nanodisc": False, "ipr": ["IPR1"], "methods": "ELECTRON MICROSCOPY",
         "title": "Cryo-EM structure of a sodium channel", "n_protein": 1, "n_unmapped": 0}
    f.update(kw)
    return f


def _ctx(**kw):
    c = SimpleNamespace(protein_dates={"P1": ["2020-01-01"]}, protein_example={"P1": "1ABC"},
                        fam_stats={"IPR1": [9, 1]}, tcdb={"P1": ["1.A.1.1.1"]})
    c.__dict__.update(kw)
    return c


def test_isotonic_is_monotone_and_inside_0_1():
    preds = [(i / 100, int(i > 60)) for i in range(100)] + [(0.95, 0), (0.05, 1)]
    knots = _isotonic(preds, bins=10)
    vals = [y for _, y in knots]
    assert vals == sorted(vals) and 0 < vals[0] and vals[-1] < 1
    assert _apply_isotonic(knots, 0.0) == vals[0] and _apply_isotonic(knots, 1.0) == vals[-1]
    assert _apply_isotonic(knots, 0.2) < _apply_isotonic(knots, 0.9)


def test_features_use_only_earlier_mpstruc_entries():
    x, why = features(_facts(), _ctx(), prior=0.1)
    assert x["same_protein"] == 1 and why["same_protein"] == "1ABC"
    x, _ = features(_facts(released="2019-06-01"), _ctx(), prior=0.1)
    assert x["same_protein"] == 0  # mpstruc only listed it later


def test_family_vote_excludes_the_structure_itself():
    x_own, _ = features(_facts(), _ctx(fam_stats={"IPR1": [1, 0]}), prior=0.1, own_label=1)
    assert x_own["family"] == 0  # its own count removed, nothing left
    x_other, _ = features(_facts(), _ctx(fam_stats={"IPR1": [1, 0]}), prior=0.1)
    assert x_other["family"] > 0


def test_tcdb_and_explanation():
    model = {"bias": -3.0, "weights": {"kw_tm": 3.0, "detergent": 1.0, "tcdb_transporter": 1.0, "same_protein": 1.0,
                                       "kw_KW-0964": -2.0}}
    x, why = features(_facts(), _ctx(), prior=0.1)
    assert x["tcdb_transporter"] == 1
    p = predict(x, model)
    reasons = explain(x, why, _facts(), p, model)
    assert p > 0.5 and reasons[0].startswith("UniProt lists")
    assert any("detergent (LMN)" in r for r in reasons) or any("TCDB" in r for r in reasons)
    # A secreted protein with no membrane keyword reads the other way round.
    f = _facts(tm_uniprots=[], keywords=["KW-0964"], detergents=[])
    x, why = features(f, _ctx(tcdb={}), prior=0.1)
    p = predict(x, model)
    reasons = explain(x, why, f, p, model)
    assert p < 0.5 and any("crossing the membrane" in r for r in reasons)
