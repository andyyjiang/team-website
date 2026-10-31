# Automation backtest

The automation is run as if it were the start of each period, knowing only the mpstruc entries released before then, and its proposals are compared with what the curators actually did. Rules, confidences, the likelihood model and the default cut-off are fitted on the tuning period only.

## Tuning period: structures released 2025-01-01 to 2025-06-30

- Released in the PDB: 8,622
- Of these, later added to mpstruc: 367
- **Flagged by the automation: 366 (99.7%)**
- Flagged but not (yet) in mpstruc: 2,418

Ranking by likelihood (how many a curator would look at each week):

| Show if likelihood ≥ | Per week | Finds (recall) | Of those shown, in mpstruc |
|---|---|---|---|
| 0.0 | 108.3 | 99.7% | 13.1% |
| 0.1 | 74.6 | 99.7% | 19.1% |
| 0.2 | 58.9 | 98.9% | 24.0% |
| 0.3 | 48.0 | 97.3% | 28.9% |
| 0.3 ← default | 46.3 | 96.2% | 29.6% |
| 0.4 | 40.1 | 92.9% | 33.1% |
| 0.5 | 32.4 | 86.4% | 38.0% |
| 0.6 | 24.9 | 77.4% | 44.3% |
| 0.7 | 15.9 | 61.9% | 55.6% |
| 0.8 | 9.9 | 43.6% | 63.0% |
| 0.9 | 3.1 | 18.0% | 83.5% |

- **Placement right (all): 80.3%** of 366
- **Placement right when the automation decides: 85.7%** of 343; 6.3% sent to curators with no default
- Subfamily right (first choice): 85.3%; in top 3: 94.3%

| Rule | Structures | Type right | Subfamily right | Parent/owner right |
|---|---|---|---|---|
| new-protein | 34 | 64.7% | 73.5% | — |
| no-evidence | 3 | 0.0% | None% | — |
| partial-protein | 20 | 0.0% | 65.0% | 25.0% |
| same-paper | 5 | 80.0% | 100.0% | 80.0% |
| same-paper-batch | 250 | 89.2% | 85.0% | 79.0% |
| same-protein | 54 | 83.3% | 100.0% | 81.5% |

Confusion (true → predicted): daughter→daughter: 45, daughter→parent: 3, daughter→related: 16, daughter→review: 8, parent→parent: 22, parent→related: 12, parent→review: 11, related→daughter: 9, related→parent: 9, related→related: 227, related→review: 4

## Test period (untouched by tuning): structures released 2025-07-01 to 2025-12-31

- Released in the PDB: 8,978
- Of these, later added to mpstruc: 447
- **Flagged by the automation: 443 (99.1%)**
- Flagged but not (yet) in mpstruc: 2,490

Ranking by likelihood (how many a curator would look at each week):

| Show if likelihood ≥ | Per week | Finds (recall) | Of those shown, in mpstruc |
|---|---|---|---|
| 0.0 | 112.2 | 99.1% | 15.1% |
| 0.1 | 80.4 | 98.9% | 21.0% |
| 0.2 | 60.6 | 96.2% | 27.2% |
| 0.3 | 49.5 | 92.8% | 32.1% |
| 0.3 ← default | 47.1 | 91.1% | 33.1% |
| 0.4 | 34.8 | 80.8% | 39.7% |
| 0.5 | 21.0 | 52.8% | 43.0% |
| 0.6 | 13.0 | 34.0% | 44.6% |
| 0.7 | 6.2 | 17.9% | 49.4% |
| 0.8 | 2.0 | 5.8% | 50.0% |
| 0.9 | 0.4 | 1.8% | 80.0% |

- **Placement right (all): 70.7%** of 443
- **Placement right when the automation decides: 79.6%** of 393; 11.3% sent to curators with no default
- Subfamily right (first choice): 85.3%; in top 3: 95.7%

| Rule | Structures | Type right | Subfamily right | Parent/owner right |
|---|---|---|---|---|
| new-protein | 42 | 69.0% | 76.2% | — |
| no-evidence | 2 | 0.0% | None% | — |
| partial-protein | 48 | 0.0% | 70.8% | 42.9% |
| same-paper | 14 | 71.4% | 100.0% | 50.0% |
| same-paper-batch | 279 | 83.9% | 85.3% | 70.5% |
| same-protein | 58 | 69.0% | 100.0% | 65.5% |

Confusion (true → predicted): daughter→daughter: 40, daughter→parent: 3, daughter→related: 32, daughter→review: 28, parent→daughter: 3, parent→parent: 29, parent→related: 17, parent→review: 8, related→daughter: 15, related→parent: 10, related→related: 244, related→review: 14

## Likelihood model (fitted on the tuning period)

| Signal | Weight |
|---|---|
| UniProt: transmembrane protein | +1.12 |
| Same paper as another candidate this week | +0.94 |
| In PDBTM (transmembrane proteins) | +0.91 |
| In MemProtMD | +0.71 |
| Same paper as an mpstruc entry | +0.64 |
| Same protein as an mpstruc entry | +0.39 |
| UniProt: membrane protein | +0.34 |
| Cryo-EM | +0.33 |
| Protein family seen in mpstruc | +0.12 |
| NMR | +0.03 |
| X-ray | -0.80 |

Title words that raise the likelihood most: agonist (+1.26), inward-facing (+1.21), methyl (+1.08), coli (+0.98), inhibitor (+0.79), lipid (+0.78), scfv16 (+0.75), lmng (+0.74), human (+0.68), obtained (+0.63), class (+0.61), msba (+0.53), xpr1 (+0.50), protein (+0.50), complexed (+0.48)

Title words that lower it most: sars-cov-2 (-1.48), map (-1.05), assembly (-1.00), spike (-0.97), nanodisc (-0.83), virus (-0.82), compound (-0.81), ca2 (-0.81), dimer (-0.77), variant (-0.61), rbd (-0.59), glycoprotein (-0.57), protease (-0.56), supercomplex (-0.54), antibody (-0.54)

## Confidence used for each rule (from the tuning period)

- parent: 48%
- partial: 0%
- related: 80%
- related-batch: 76%
- daughter: 83%
