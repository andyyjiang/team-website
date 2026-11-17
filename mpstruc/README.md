# mpstruc curation tool

A database browser and automation pipeline for
[mpstruc](https://blanco.biomol.uci.edu/mpstruc/), *Membrane Proteins of Known 3D Structure*
(created at UC Irvine by Stephen White and Craig Snider, and now jointly managed by members
of the USC Bridge Institute). This is the curation team's tool.

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
python3 -m mpstruc.automation apply-rules  # curator_rules.json on the current candidates, no weekly run
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
| `curator_rules.json` | The curators' own placement rules, one item per subfamily (see "Curators' rules" below) |
| `mpstruc/rules.py` | Loads and checks `curator_rules.json`; matches each rule to a subfamily by name |
| `mpstruc/barrel.py` | Counts the strands of a membrane beta-barrel from a structure's C-alpha coordinates |
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

## Curators' rules

Curators can state how a subfamily is recognised, and the automation follows it. The rules are
in `curator_rules.json`; add one by adding an item and running
`python3 -m mpstruc.automation apply-rules` (the weekly job reads the file too):

```json
{"subgroup": "Omp85-TpsB Outer Membrane Transporter Superfamily", "subgroupKey": "G02-S007",
 "note": "Typically 16-stranded beta-barrel.", "barrelStrands": 16, "tolerance": 1,
 "addedBy": "curator", "addedOn": "2026-11-13"}
```

- `subgroup` is the subfamily's name exactly as in mpstruc (rules are matched by name; the key
  is kept for reference). `note` is the curator's sentence, shown on the subfamily's page.
- `barrelStrands` (optional) is the typical number of strands of the membrane barrel, and
  `tolerance` how far a counted number may be from it. Other criteria can be added as new fields.
- **Counting.** `mpstruc/barrel.py` counts the strands from the C-alpha coordinates of the
  outer-membrane protein's chains (PDB ModelServer, cached): strands that are hydrogen-bonded side
  by side in a ring, in chain order. POTRA, plug and passenger domains are other sheets and are
  not counted; barrels made by several chains (trimeric autotransporters, MspA) are counted over the
  assembly. The count of a barrel whose seam is open (BamA or Sam50 with a substrate) is a lower
  bound and is shown as "at least N"; when two rings fit the same strands it is shown as a range
  ("8 to 10"). `python3 -m mpstruc.barrel validate` compares it with known structures
  (`reports/barrel_validation.md`: 56 of 57 exactly right, 14 of 16 open barrels).
- **Membrane barrels only.** UniProt's keywords describe the full-length protein, so a ring of
  strands is read as a membrane barrel only if it lies in the membrane part UniProt names
  (transmembrane strands, a barrel domain, the membrane anchor or translocator domain); a ring
  made by several chains that cannot be checked this way needs 12 or more strands. The head and
  stalk fragments of trimeric autotransporters, which close rings of 6 to 17 strands, are in the
  guard set of `reports/barrel_validation.md`, and the tests fail if one is read as a barrel.
- **Placement.** If the count contradicts the suggested subfamily's rule and fits another ruled
  subfamily, the suggestion moves there (TamA, 16 strands: Omp85-TpsB, not Autotransporters). The
  barrel that is checked is the best one known for the paper: a closed barrel in the structure,
  else in another structure of the same paper, else UniProt's membrane strands if all are
  resolved; and a structure that would join an official entry which itself contradicts the rule
  (TamA 4C00) moves with it. An open barrel decides only what its lower bound decides: 13 strands
  seen in an open barrel are not "a 12-strand autotransporter". A barrel made by several chains
  moves nothing. A new protein that only its protein families put in a ruled subfamily, counted
  without a barrel and with a membrane score under 50%, goes to review; a structure whose
  coordinates could not be fetched is left as it is ("rule not checked"), never treated as
  having no barrel. Each candidate records `barrelStrands` and `ruleCheck`.
- **Membrane score.** A barrel of 8 or more strands counted in an outer-membrane protein counts as
  a resolved barrel when UniProt names no barrel domain (BamA in BAM complexes). The model is not
  retrained; `python3 -m mpstruc.membrane barrels` measures the effect (`reports/barrel_membrane.md`).
- **Official entries.** `reports/rule_conflicts.md` and `../data/curator_rules.json` list the
  official mpstruc entries whose barrel contradicts a rule. This is information for the curators;
  the automation never changes official mpstruc.

## How the automation decides

Rules come from `reports/placement_analysis.md`. They were tested by replaying July to December 2025
week by week, with each structure as it looked in its release week (no paper linked yet if the PDB
linked it later, no PDB annotations added later), against what curators did (`reports/backtest.md`):

| A new structure that… | is suggested as | fully right* |
|---|---|---|
| comes from the same paper as an mpstruc entry, same protein | related entry of that entry | 67% |
| comes from the same paper as other new structures, same protein | related to the lowest PDB ID among them | 66% |
| is the same protein(s) as an mpstruc entry, new paper | daughter of the parent holding most of that protein | 69% |
| shares only some proteins (complex or fusion differs) | curator decides (no default) | — |
| is a protein mpstruc does not have | new parent in the best-matching subfamily | 46% |

\* type, subfamily and parent/owner all match the curators' choice. Overall, at the default cut-off,
the automation finds 92% of what mpstruc adds, shows about 50 structures a week (about 1 in 3 of
them are added), picks the right subfamily first 93% of the time, and gets the whole placement right
64% of the time (type alone: 81%).

- **Membrane proteins in a structure** are the UniProt accessions with a transmembrane keyword, not
  the PDB's OPM/PDBTM/MemProtMD tags, which new releases no longer receive.
- **"Same paper"** uses the PDB's citation, or, until a paper is linked (about half of new releases
  in their first weeks), the first and last depositing author.
- **Families** come from UniProt's InterPro entries, which exist on release day.
- **The likelihood that mpstruc lists a structure** comes from a small logistic model over these
  signals, the experimental method and words in the PDB title. Its cut-off is chosen by
  cross-validation grouped by paper. RCSB's own copy of the mpstruc classification is never used.
- **Limitations:**
  - Curation lag is not modelled; shifting the cut-off date by 30–90 days moved the results by
    1 point or less.
  - UniProt families and keywords are today's versions.
  - The rules were designed on the full snapshot, which includes the test period.
- **Curator decisions** are saved in the browser and leave it only as downloaded CSV or mpstruc XML
  blocks. The weekly job re-checks the last 26 weeks, so suggestions improve as papers are linked.
  Every other PDB entry citing a candidate's paper, released before or after those 26 weeks and not in
  mpstruc or decided, is taken in too (`mpstruc/siblings.py`, marked `sibling` in `data/candidates.json`) when
  it is at least 10% likely to be a membrane protein, or is the same protein as a kept candidate of the paper,
  so a study is suggested whole under its first-deposited structure (reissued IDs such as 24UG count as newer
  than every 9xxx ID); `python -m mpstruc.automation siblings` adds them to the current file without a weekly
  run. The paper lookup is cached for 3 days; if the PDB search fails, a paper keeps its last cached entries,
  or in CI (no cache) the entries the previous `data/candidates.json` held.
