import json

import pytest

from mpstruc import markup
from mpstruc.markup import LEXICON_PATH, format_text, strip_subsup


@pytest.mark.parametrize("plain, marked", [
    ("Ca2+", "Ca<sup>2+</sup>"),
    ("Na+/K+-ATPase", "Na<sup>+</sup>/K<sup>+</sup>-ATPase"),
    ("cytochrome b6f complex", "cytochrome b<sub>6</sub>f complex"),
    ("cytochrome bc1 complex", "cytochrome bc<sub>1</sub> complex"),
    ("F1Fo ATP synthase", "F<sub>1</sub>F<sub>o</sub> ATP synthase"),
    ("FoF1", "F<sub>o</sub>F<sub>1</sub>"),
    ("GABAA receptor", "GABA<sub>A</sub> receptor"),
    ("5-HT2A serotonin receptor - miniGq complex",
     "5-HT<sub>2A</sub> serotonin receptor - miniG<sub>q</sub> complex"),
    ("Glucagon receptor-Gs complex", "Glucagon receptor-G<sub>s</sub> complex"),
    ("Super Complex III2IV", "Super Complex III<sub>2</sub>IV"),
    ("CIII2", "CIII<sub>2</sub>"),
    ("PI(3,5)P2", "PI(3,5)P<sub>2</sub>"),
    ("H2O", "H<sub>2</sub>O"),
    ("10 mM MgCl2", "10 mM MgCl<sub>2</sub>"),
    ("NH4+", "NH<sub>4</sub><sup>+</sup>"),
    ("HCO3- transporter", "HCO<sub>3</sub><sup>-</sup> transporter"),
    ("SO42-", "SO<sub>4</sub><sup>2-</sup>"),
    ("Cl- channel", "Cl<sup>-</sup> channel"),
    ("Sr2+-bound", "Sr<sup>2+</sup>-bound"),          # general cation rule
    ("with bound ADP.BeF3", "with bound ADP.BeF<sub>3</sub>"),
    ("AlF4-", "AlF<sub>4</sub><sup>-</sup>"),
    ("KH2PO4 buffer", "KH<sub>2</sub>PO<sub>4</sub> buffer"),
    ("Mn4CaO5 cluster", "Mn<sub>4</sub>CaO<sub>5</sub> cluster"),
    # Two ions: not PDB IDs, even though they start with a digit.
    ("SERCA in E1-2Ca2+ state", "SERCA in E1-2Ca<sup>2+</sup> state"),
    ("2NO3− symport", "2NO<sub>3</sub><sup>−</sup> symport"),
    # Short keys in the context curators use them in.
    ("Histamine H4 receptor", "Histamine H<sub>4</sub> receptor"),
    ("leukotriene B4 receptor", "leukotriene B<sub>4</sub> receptor"),
    ("Neuropeptide Y Y1 receptor", "Neuropeptide Y Y<sub>1</sub> receptor"),
    ("endothelin ETA receptor", "endothelin ET<sub>A</sub> receptor"),
    ("space group P63", "space group P6<sub>3</sub>"),
    # G-protein pairs: mpstruc subscripts the whole pair (miniG<sub>s/q</sub>).
    ("miniGs/q complex", "miniG<sub>s/q</sub> complex"),
    ("Gi/o and Gq/11", "G<sub>i/o</sub> and G<sub>q/11</sub>"),
    ("G12/13 coupling", "G<sub>12/13</sub> coupling"),
    # "Go" at the start is the protein when a hyphen follows.
    ("Go-Gi chimera", "G<sub>o</sub>-G<sub>i</sub> chimera"),
])
def test_marks_up_like_mpstruc(plain, marked):
    assert format_text(plain) == marked


@pytest.mark.parametrize("text", [
    "PI(4,5)P2",           # mpstruc mostly leaves this plain
    "PIP2", "IP3", "PGH2", "KCNK3", "FKS1", "IF1", "CNIH2",
    "Cl-dependent",        # a hyphen, not a charge
    "Done. Go home",       # English word at sentence start
    "PDB 1H2O and 9ABC",   # PDB IDs
    "see https://example.org/Ca2+/H2O",
    "3.10 &Aring;: 8TSI",
    "Calcium ATPase",      # no tag inside words
    # Gene/protein names and ligand codes made of element letters.
    "PHO84 phosphate transporter", "OPN4", "OSH6", "Protein SOP4", "SCO1",
    "NOS1", "COS7 cells", "CHO1 synthase", "CO1", "POS5", "NaCl1",
    "the IST2-OSH6 complex", "with bound CO89",
    # Short keys outside the context they were learned from.
    "TYR-DAL-GLY-MEA-ETA", "Histone H4 tail", "N-terminal H1 helix",
    "Y1 phosphorylation", "chain P63 residue", "helix H1 and H4", "subunit B4 of",
    "",
])
def test_leaves_alone(text):
    assert format_text(text) == text


def test_none_passthrough():
    assert format_text(None) is None


def test_existing_markup_untouched_and_idempotent():
    s = "Ca<sup>2+</sup>-bound, Ca2+ free, G<sub>s</sub> and Gs"
    once = format_text(s)
    assert once == "Ca<sup>2+</sup>-bound, Ca<sup>2+</sup> free, G<sub>s</sub> and G<sub>s</sub>"
    assert format_text(once) == once


def test_partial_markup_not_double_tagged():
    # A word already partly marked up is left exactly as written.
    for s in ("F<sub>1</sub>Fo", "Ca2<sup>+</sup>", "b<sub>6</sub>f", "GABA<sub>A</sub>"):
        assert format_text(s) == s
    # Tag attributes are never edited, the link text is.
    assert format_text("<a href='x?Ca2+'>Ca2+</a>") == "<a href='x?Ca2+'>Ca<sup>2+</sup></a>"


def test_br_is_a_boundary():
    assert format_text("cryo-EM structure<br /> with Ca2+, 3.1 Å: 9ABC") == \
        "cryo-EM structure<br /> with Ca<sup>2+</sup>, 3.1 Å: 9ABC"


def test_greek_entity_and_character_both_work():
    lex = {"α2A": "α<sub>2A</sub>"}
    assert format_text("α2A adrenergic", lex) == "α<sub>2A</sub> adrenergic"
    assert format_text("&alpha;2A adrenergic", lex) == "&alpha;<sub>2A</sub> adrenergic"


def test_lexicon_file_is_current_and_curated():
    data = json.loads(LEXICON_PATH.read_text(encoding="utf-8"))
    lex = data["lexicon"]
    assert lex["Ca2+"] == "Ca<sup>2+</sup>" and lex["Gs"] == "G<sub>s</sub>"
    for plain, marked in lex.items():
        assert strip_subsup(marked) == plain
        st = data["stats"][plain]
        assert st["consistency"] >= data["min_consistency"] and st["kept"]
    # Spellings mpstruc curators write both ways are left out.
    for ambiguous in ("PIP2", "IP3", "PI(4,5)P2", "β2"):
        assert ambiguous not in lex


def test_units_mining():
    got = list(markup.units("Na<sup>+</sup>/K<sup>+</sup>-ATPase and cytochrome b<sub>6</sub>f "
                            "and &beta;<sub>2</sub> and PI(3,5)P<sub>2</sub>"))
    assert ("Na+", "Na<sup>+</sup>") in got
    assert ("b6f", "b<sub>6</sub>f") in got
    assert ("β2", "β<sub>2</sub>") in got
    assert ("PI(3,5)P2", "PI(3,5)P<sub>2</sub>") in got


def test_whole_corpus_safe():
    """On mpstruc's own fields: idempotent, existing text never altered, and
    fields without any markup are almost always left unchanged."""
    recs = list(markup.corpus([markup.FIXTURES / markup.SNAPSHOTS[1]]))
    changed_plain = plain = 0
    for _, field, raw in recs:
        out = format_text(raw)
        assert strip_subsup(out) == strip_subsup(raw)
        assert format_text(out) == out
        if field in markup.CURATOR_FIELDS and "<su" not in raw:
            plain += 1
            changed_plain += out != raw
    assert changed_plain / plain < 0.01


def test_held_out_new_entries():
    """The 19 entries curators added between the two snapshots, with the
    lexicon learned only from the older snapshot."""
    old = list(markup.corpus([markup.FIXTURES / markup.SNAPSHOTS[0]]))
    lex = markup.build_lexicon(old)
    new = [r for r in markup.corpus([markup.FIXTURES / markup.SNAPSHOTS[1]])
           if r[0] in markup.HELD_OUT_19 and r[1] in markup.CURATOR_FIELDS]
    res = markup.score(new, lex)
    assert res["recall"] == 1.0
    assert res["correct"] >= res["tags_added"] - 1   # only 8FYN "A2A" (curator left it plain)
