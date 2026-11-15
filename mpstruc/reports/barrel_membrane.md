# Membrane score: barrels counted from coordinates

A barrel of at least 8 strands counted in a chain of an outer-membrane protein (UniProt: transmembrane
beta strand, cell outer membrane and not only a lipoprotein, or an Omp85 family) now counts as a resolved
membrane barrel when UniProt maps neither strands nor a barrel domain for the protein, and the ring lies in
the protein's membrane part (`membrane.membrane_barrel`). The model was not retrained. Measured on the test
period (2024-09-01 to 2025-03-31) with `python -m mpstruc.membrane barrels`.

- Structures in the period: 9583; with such a protein: 16; the counted barrel raised the score of 8 (it is not used where it would lower one).
- Crossed 50%: 2 (2 labelled membrane proteins, 0 labelled not, 0 left out of the labels).
- **What this period can and cannot show.** Only 16 of its structures have a protein the count applies to, so "no false positive here" says little about soluble fragments of outer-membrane
  proteins. Those are tested on purpose: the guard set in `reports/barrel_validation.md` (heads and stalks of
  trimeric autotransporters, POTRA-only, passenger, usher and TonB-dependent receptor soluble domains, lipoprotein
  barrels) must give no membrane barrel, and `tests/test_rules.py` fails if the report says otherwise.
- Structures below that are labelled "not membrane" although a 22-strand TonB-dependent receptor barrel is in them
  (FoxA, FemA) are wrong in the labels (PDBTM and MemProtMD had not listed them yet), not in the score; they
  lower the precision reported here before and after alike.

| | Called membrane (>= 50%) | Precision | Recall | AUC | Brier |
|---|---|---|---|---|---|
| Before | 1364 | 84.0% | 96.5% | 0.9855 | 0.0235 |
| After | 1366 | 84.0% | 96.7% | 0.9854 | 0.0236 |

## Recall by kind (share scored at least 50%)

| Kind | n | Before | After |
|---|---|---|---|
| MONOTOPIC MEMBRANE PROTEINS | 7 | 0.0% | 0.0% |
| TRANSMEMBRANE PROTEINS: ALPHA-HELICAL | 474 | 96.4% | 96.4% |
| TRANSMEMBRANE PROTEINS: BETA-BARREL | 19 | 89.5% | 89.5% |
| outer membrane (UniProt) | 32 | 84.4% | 90.6% |

## Structures in which a barrel was counted

| PDB | Strands | Before | After | Label | In mpstruc | Title |
|---|---|---|---|---|---|---|
| 8RM3 | 22 | 59% | 92% | not membrane | no | Crystal structure of ferrioxamine transporter FoxA |
| 8RMI | 22 | 59% | 92% | not membrane | no | Crystal structure of ferrioxamine transporter FoxA |
| 8S34 | 22 | 64% | 93% | not membrane | no | Ferric-mycobactin receptor (FemA) in complex with aeruginic acid |
| 8XNB | 16 | 16% | 70% | membrane | no | Cryo-EM structure of Deinococcus radiodurans BamA |
| 9BOG | 24 | 24% | 73% | membrane | no | Structural basis for adhesin secretion by the outer-membrane usher in type 1 pili |
| 9EX3 | 22 | 64% | 93% | membrane | no | Ferric-mycobactin receptor (FemA) in complex with dihydroaeruginoic acid |
| 9F2T | 22 | 64% | 93% | membrane | no | Ferric-mycobactin receptor (FemA) |
| 9FVQ | 22 | 64% | 93% | not membrane | no | Ferric-mycobactin receptor (FemA) in complex with pyochelin |
