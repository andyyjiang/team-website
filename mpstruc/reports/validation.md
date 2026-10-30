# Deliverable 1: parser validation

## Summary

- **Round trip:** parsing the fixture and writing it back is **byte-identical**.
- **Totals:** all match the xlsx `Source` sheet (table below).
- **Cell-by-cell diff against the xlsx:** 0 unexplained difference(s). 15 explained ones are listed below.
- **xlsx defect:** the `Structures` sheet has 16 blank rows; the 16 records missing from it are the whole G01-S001-P0001 branch (1PRH, ram COX-1). They are present in the XML.
- **Integrity findings:** 10 errors, 209 warnings, 147 info (see `integrity.md`).

## Snapshot identity

- Fixture SHA-256: `4174f67bb5685afe0c0cf046f0c5bc5ceaa53162f8fd9568cf679f5187ea9f7a`
- SHA-256 recorded in the xlsx `Source` sheet: `1ce45f6a3176467bfec76000d594efa3a88196d0169c43460253a7e5cce2e58b`. These differ. The timestamps match and every count below matches, so the xlsx hash was most likely taken over different bytes (for example before re-saving). The fixture hash above is the one the parser checks.

## Totals

| Quantity | Parsed from XML | xlsx `Source` sheet | Match |
|---|---|---|---|
| Parent records | 2,070 | 2,070 | ✓ |
| Daughter records | 2,248 | 2,248 | ✓ |
| Related PDB listings | 6,273 | 6,273 | ✓ |
| All PDB listings | 10,591 | 10,591 | ✓ |
| Subfamilies | 206 | 206 | ✓ |
| Distinct parent PDB IDs | 2,064 | 2,064 | ✓ |
| Distinct PDB IDs in all listings | 10,552 | 10,552 | ✓ |
| Bibliography records | 4,331 | 4,331 | ✓ |

## Structures sheet

- Parsed listings: 10,591; populated xlsx rows: 10,575; blank xlsx rows: 16
- Rows in the xlsx that we do not produce: 0
- Records missing from the xlsx: 16

Missing records, grouped by parent:

- `G01-S001-P0001`: `G01-S001-P0001`, `G01-S001-P0001-D001`, `G01-S001-P0001-D002`, `G01-S001-P0001-D003`, `G01-S001-P0001-D003-R001`, `G01-S001-P0001-D003-R002`, `G01-S001-P0001-D003-R003`, `G01-S001-P0001-D004`, `G01-S001-P0001-D005`, `G01-S001-P0001-D006`, `G01-S001-P0001-D007`, `G01-S001-P0001-D007-R001`, `G01-S001-P0001-D007-R002`, `G01-S001-P0001-D007-R003`, `G01-S001-P0001-D007-R004`, `G01-S001-P0001-D008`

### Explained cell differences

| Class | Meaning | Count |
|---|---|---|
| resolution_parsed | The xlsx left a parseable resolution blank (e.g. `1.77 Å`, `1.8 Å (P63 crystal form)`, `3.4/3.7 Å`) | 7 |
| entity_repaired | The xlsx kept a misspelled entity (`&Ariing;`, `&alfa;`, `& Aring;`, `Aring;`) as literal text; we display Å / α | 6 |
| invisible_removed | The xlsx kept zero-width characters | 2 |

| Record key | Column | Parsed | xlsx | Class |
|---|---|---|---|---|
| `G02-S003-P0029` | Description | NHis-native, 1.83 Å: 4V3H ⏎ wild-type, 2.30 Å: 4D51 ⏎ α-cyclodextrin soak, 1.70 Å: 4D5B ⏎ β-cyclodextrin soak, 1.90 Å: 4D5D | NHis-native, 1.83 Å: 4V3H ⏎ wild-type, 2.30 Å: 4D51 ⏎ &alfa;-cyclodextrin soak, 1.70 Å: 4D5B ⏎ β-cyclodextrin soak, 1.90 Å: 4D5D | entity_repaired |
| `G02-S005-P0054-D003` | Resolution (Å) | 1.77 |  | resolution_parsed |
| `G03-S033-P0001-D001` | Description | cryo-EM structure ⏎ Ca+2-bound intermediate state, 3.7 Å: 6QMA ⏎ Ca+2-bound closed state, 3.6 Å:6QMB ⏎ Ca+2-free state, 3.8 Å: 6QM4 ⏎ in … | cryo-EM structure ⏎ Ca+2-bound intermediate state, 3.7 Å: 6QMA ⏎ Ca+2-bound closed state, 3.6 & Aring;:6QMB ⏎ Ca+2-free state, 3.8 Å: 6QM… | entity_repaired |
| `G03-S037-P0010-D001` | Description | x-ray structure ⏎ in complex with an inhibitory Nb61, 3.50 Å: 8QZ2 ⏎ in complex with an activatory Nb67, 2.40 Å: 8QZ3 ⏎ in complex with a… | x-ray structure ⏎ ​​in complex with an inhibitory Nb61, 3.50 Å: 8QZ2 ⏎ in complex with an activatory Nb67, 2.40 Å: 8QZ3 ⏎ in complex with… | invisible_removed |
| `G03-S038-P0001-D001` | Description | [Ca2+] = 5mM ⏎ W195Y mutant in complex with Br-dihydropyridine derivative UK-59811, 3.3 Å: 5KLG ⏎ wild-type protein in complex with Br-di… | [Ca2+] = 5mM ⏎ W195Y mutant in complex with Br-dihydropyridine derivative UK-59811, 3.3 Å: 5KLG ⏎ wild-type protein in complex with Br-di… | entity_repaired |
| `G03-S047-P0001-D001` | Resolution (Å) | 1.8 |  | resolution_parsed |
| `G03-S051-P0035` | Resolution (Å) | 3.4 |  | resolution_parsed |
| `G03-S051-P0165` | Resolution (Å) | 2.91 |  | resolution_parsed |
| `G03-S054-P0014` | Resolution (Å) | 3.4 |  | resolution_parsed |
| `G03-S070-P0003` | Resolution (Å) | 2.8 |  | resolution_parsed |
| `G03-S094-P0001-D011` | Description | This and the additional structures below reveal two discrete multisite binding pockets. ⏎ Unliganded AcrB, 3.35 Å: 3AOA ⏎ With bound eryt… | This and the additional structures below reveal two discrete multisite binding pockets. ⏎ Unliganded AcrB, 3.35 Å: 3AOA ⏎ With bound eryt… | entity_repaired |
| `G03-S133-P0021-D002` | Description | cryo-EM structure ⏎ in the E2P transition state, 3.98 Å: 7KY5 ⏎ in the E2P state, 3.50 Å: 7KYA ⏎ Dnf1-Lem3 complex in the apo E1 state, 3… | cryo-EM structure ⏎ in the E2P transition state, 3.98 Å: 7KY5 ⏎ in the E2P state, 3.50 &Ariing;: 7KYA ⏎ Dnf1-Lem3 complex in the apo E1 s… | entity_repaired |
| `G03-S140-P0019` | Description | cryo-EM structure ⏎ in the dimer-of-dimer state, 2.70 Å: 6WXU ⏎ complex in the presence of NADPH, 3.30 Å: 6WXV | cryo-EM structure ⏎ in the dimer-of-dimer state, 2.70 Å: 6WXU ⏎ complex in the presence of NADPH, 3.30 &Ariing;: 6WXV | entity_repaired |
| `G03-S143-P0003` | Resolution (Å) | 1.6 |  | resolution_parsed |
| `G03-S157-P0002-D001` | Description | x-ray structure ⏎ soaked with CS2, 2.65 Å: 8STD | x-ray structure ⏎ ​​soaked with CS2, 2.65 Å: 8STD | invisible_removed |

## References sheet

- Parsed bibliography records: 4,331; xlsx rows: 4,331
- Missing / extra / differing: 0 / 0 / 0
- The 9 reference rows for the 1PRH branch are present in this sheet even though the branch's structure rows are missing from `Structures`.

## Subfamilies sheet

- Parsed subgroups: 206; xlsx rows: 206; differing cells: 0
- Subgroups are numbered within their group (`G03-S136`), matching the xlsx record keys. Two subgroups are named *Hydrolases* (G01-S005, G03-S136), so subgroups are identified by key, never by name.

## Known limitations of the xlsx (kept in our data, lost in the xlsx)

- Formatting is stripped: `P4<sub>3</sub>2<sub>1</sub>2` becomes `P43212`, and `Mn<sup>2+</sup>` becomes `Mn2+`. We keep the raw markup and render it.
- Related entries carry no metadata in the XML. We will fill it in from the PDB at display time, labelled "from PDB".
