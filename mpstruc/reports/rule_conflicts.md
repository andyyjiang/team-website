# Official mpstruc entries that contradict a curators' rule

Barrel strands counted from the coordinates (`mpstruc/barrel.py`) against the rules in `curator_rules.json`.
For the curators' information: the automation cannot change official mpstruc. Regenerate with
`python -m mpstruc.automation apply-rules`.

## Rules

- **Outer Membrane Autotransporters** (G02-S006): Typically 12-stranded beta-barrel. (added by curator, 2026-11-13; counts within 1 of 12 fit)
- **Omp85-TpsB Outer Membrane Transporter Superfamily** (G02-S007): Typically 16-stranded beta-barrel. (added by curator, 2026-11-13; counts within 1 of 16 fit)

## Entries that do not fit (2)

| PDB | Entry | Listed under | Rule | Counted | Fits instead |
|---|---|---|---|---|---|
| 6WIL | CdiB autotransporter (G02-S006-P0013) | Outer Membrane Autotransporters | 12 strands | 16 strands (chain A) | Omp85-TpsB Outer Membrane Transporter Superfamily |
| 6WIM | CdiB autotransporter (G02-S006-P0014) | Outer Membrane Autotransporters | 12 strands | 16 strands (chain A) | Omp85-TpsB Outer Membrane Transporter Superfamily |
