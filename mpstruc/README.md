# mpstruc curation tool

A database browser and automation pipeline for
[mpstruc](https://blanco.biomol.uci.edu/mpstruc/), *Membrane Proteins of Known 3D Structure*
(created by Stephen White, maintained by Craig Snider, UC Irvine). This is a tool for the
mpstruc curators; mpstruc itself remains the authoritative source.

## Status

| Deliverable | State |
|---|---|
| 1. XML parser, validation against the reference xlsx, integrity report | **done, awaiting review** |
| 2. Database tab (`../mpstruc.html`) | done |
| 3. Placement analysis, automation pipeline, evaluation harness | done (`reports/placement_analysis.md`, `reports/backtest.md`) |
| 4. Automation tab (`../mpstruc.html#/automation`) and weekly job | done, in review |
| 5. Final README | this file grows with each deliverable |

## Quick start

```bash
pip install openpyxl pytest
python3 -m pytest -q              # 46 tests, offline
python3 -m mpstruc.report         # writes reports/ (uses RCSB; add --offline to skip)
python3 -m mpstruc.web            # writes ../data/mpstruc.json for the website (uses RCSB + PubMed)
python3 -m mpstruc.placement      # how curators place structures -> reports/placement_analysis.md
python3 -m mpstruc.automation backtest   # fit + test the automation -> model.json, reports/backtest.md
python3 -m mpstruc.automation run        # this week's suggestions -> ../data/candidates.json
```

Set `MPSTRUC_XML=/path/to/newer.xml` to use a newer mpstruc export than the saved snapshot.

The database page is `mpstruc.html` at the site root. It loads `data/mpstruc.json`; regenerate
that file with `python3 -m mpstruc.web` whenever the XML snapshot changes. Preview it with
`python3 -m http.server` in the site root and open http://localhost:8000/mpstruc.html.

## Layout

| Path | What it is |
|---|---|
| `fixtures/mpstruc_2026-10-02.xml` | The 2 Oct 2026 mpstruc XML snapshot (SHA-256 `4174f67b…`) |
| `fixtures/mpstruc_XML.xlsx` | Reference flattening of the same snapshot |
| `mpstruc/xmlparse.py` | Lossless line-based reader/writer (the export is not well-formed XML) |
| `mpstruc/model.py` | Typed model (groups → subgroups → parents → daughters → related), record keys, exact XML writer |
| `mpstruc/text.py` | Raw markup → display text, entity repair, inline PDB citations |
| `mpstruc/resolution.py` | Resolution grammar |
| `mpstruc/flatten.py` | 4-sheet xlsx layout (Structures, References, Subfamilies) |
| `mpstruc/validate.py` | Cell-by-cell comparison with the reference xlsx |
| `mpstruc/integrity.py` | Data-integrity checks (offline and RCSB-backed) |
| `mpstruc/rcsb.py` | Cached RCSB client (holdings, entry summaries) |
| `mpstruc/report.py` | Writes `reports/validation.md`, `integrity.md`, and the JSON versions |
| `mpstruc/web.py` | Exports `../data/mpstruc.json` for the website (sanitised HTML, PDB details, data-issue flags) |
| `mpstruc/pubmed.py` | Cached PubMed client used by the citation checks |
| `mpstruc/placement.py` | Measures how curators split parent / daughter / related (UniProt, organism, paper) |
| `mpstruc/automation.py` | Finds new PDB releases, scores how likely mpstruc lists them, proposes placements; backtest |
| `model.json` | Fitted likelihood model, default cut-off and per-rule confidence (from `backtest`) |
| `../.github/workflows/mpstruc-weekly.yml` | Thursday job: refresh `../data/candidates.json` and `../data/mpstruc.json` |

## Design decisions

- **Why not a standard XML parser.** The export contains HTML entities that XML doesn't
  define (`&Aring;`, `&nbsp;`), misspelled entities (`&Ariing;`, `&alfa;`), unescaped `&`, and
  broken tags (`<ahref='…'>`, an unclosed `<a>`, `<super>`). Strict parsers reject the file, and
  "recover" modes silently lose text. The file is strictly one element per line, so we parse
  that grammar and keep every field's raw content byte-for-byte. Writing the model back
  reproduces the file exactly, and a test enforces this.
- **Schema.** The element order for every record type is enforced in `model.py`. Any
  deviation raises `MpstrucParseError`. This replaces an XSD, which could not validate a file
  that isn't well-formed XML.
- **Record keys.** Keys look like `G03-S136-P0002-D001-R003`: group, then subgroup *within its
  group*, then parent, daughter, and related entry. They match the reference xlsx. Subgroups are
  always identified by key, because two subgroups are both named *Hydrolases*.
- **Raw vs display.** Stored fields are never modified. Display text comes from `to_text`:
  `<br />` becomes a newline, other tags are dropped, entities are decoded, and the known
  misspellings are repaired.
- **Resolution.** `1.8 Å (P63 crystal form)` parses to 1.8. `3.4/3.7 Å` parses to 3.4
  (kind *multiple*). NMR and electron-crystallography text becomes null with a method hint. `-`
  and `HEK293 cells` are *unparseable* and reported. Missing resolution is null, never 0.

## Reports

- `reports/validation.md`: round trip, totals, and the xlsx diff. Every difference is
  classified, and the test suite requires zero unexplained differences.
- `reports/integrity.md`: findings for curators, grouped by check. Nothing is auto-corrected.

## How the automation decides

Measured on the snapshot (`reports/placement_analysis.md`), then tested on PDB releases from
July to December 2025 that the model never saw (`reports/backtest.md`):

| A new structure that… | is suggested as | right in testing |
|---|---|---|
| comes from the same paper as an mpstruc entry, same protein | related entry of that entry | ~71% |
| comes from the same paper as other new structures, same protein | related to the lowest PDB ID among them | ~84% |
| is the same protein(s) as an mpstruc entry, new paper | daughter of the parent holding most of that protein | ~69% |
| shares only some proteins (complex or fusion differs) | curator decides (no default) | — |
| is a protein mpstruc does not have | new parent in the best-matching subfamily | ~69% |

The likelihood that mpstruc lists a structure at all comes from a small logistic model over
membrane signals (OPM, PDBTM, MemProtMD, UniProt keywords), links to mpstruc (same protein,
family, paper), the experimental method and words in the PDB title. RCSB's own copy of the
mpstruc classification is never used. Curator decisions are saved in the browser and leave it
only as downloaded CSV or mpstruc XML blocks.
