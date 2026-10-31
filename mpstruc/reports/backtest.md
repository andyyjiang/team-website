# Automation backtest

Each week of a past period is replayed: the automation knows only the mpstruc entries released before that week and is scored against what mpstruc's curators eventually did. The likelihood model, the default cut-off and the per-rule confidences are fitted on Jan–Jun 2025 only (cut-off by cross-validation grouped by paper). The placement rules themselves were designed from the full 2 Oct 2026 snapshot (reports/placement_analysis.md), which includes the test period (about 4% of it).

## Summary: test period (Jul–Dec 2025, never used for fitting), at the default cut-off

| Data the automation has | Finds | Per week | Shown that mpstruc added | Subfamily 1st (top 3) | Type right* | Fully right* | No default |
|---|---|---|---|---|---|---|---|
| release week | 92.4% | 49.7 | 31.6% | 93.0% (95.9%) | 81.2% | 64.2% | 5.9% |
| once papers are linked | 91.7% | 50.1 | 31.1% | 92.8% (95.5%) | 82.1% | 64.5% | 5.4% |

\* of the structures the automation places itself. *Release week* is what curators see first. The weekly job re-checks the last 26 weeks every Thursday, so as papers are linked to a structure its suggestion moves toward the *once papers are linked* row. About half of new releases have no paper linked in their first weeks.

## Test period, release week

Structures released 2025-07-01 to 2025-12-31, replayed one week at a time as each looked in its release week.

- Released: 8,978; later added to mpstruc: 447; flagged at any likelihood: 442 (98.9%)

| Show if likelihood ≥ | Per week | Finds (recall) | Of those shown, in mpstruc |
|---|---|---|---|
| 0.00 | 106.1 | 98.9% | 15.9% |
| 0.10 | 82.7 | 98.9% | 20.3% |
| 0.20 | 67.3 | 98.7% | 24.9% |
| 0.30 | 57.4 | 96.9% | 28.7% |
| 0.40 | 50.3 | 92.6% | 31.3% |
| 0.41 ← default | 49.7 | 92.4% | 31.6% |
| 0.50 | 39.1 | 81.2% | 35.3% |
| 0.60 | 30.7 | 69.6% | 38.5% |
| 0.70 | 15.0 | 37.6% | 42.7% |
| 0.80 | 4.8 | 14.5% | 51.6% |
| 0.90 | 0.9 | 3.1% | 60.9% |

- Placement type right: 76.5% of all 442; 81.2% of the 416 the automation decides (5.9% go to curators with no default)
- Placement fully right (type + subfamily + parent/owner): 60.4% of all; 64.2% of those it decides
- Subfamily right first time: 93.0%; in top 3: 95.9% (of all 442; a subfamily was suggested for 99.5%)

| Rule | Structures | Type right | Fully right | Subfamily right | Parent/owner right |
|---|---|---|---|---|---|
| new-protein | 52 | 63.5% | 46.2% | 80.8% | — |
| no-evidence | 2 | 0.0% | 0.0% | None% | — |
| partial-protein | 24 | 0.0% | 0.0% | 87.5% | 52.6% |
| same-paper | 18 | 77.8% | 66.7% | 100.0% | 75.0% |
| same-paper-batch | 252 | 88.1% | 65.9% | 93.7% | 75.0% |
| same-protein | 94 | 73.4% | 69.1% | 100.0% | 71.4% |

Confusion (true → predicted): daughter→daughter: 69, daughter→parent: 3, daughter→related: 20, daughter→review: 10, parent→daughter: 3, parent→parent: 33, parent→related: 14, parent→review: 7, related→daughter: 22, related→parent: 16, related→related: 236, related→review: 9

## Test period, once papers are linked

Structures released 2025-07-01 to 2025-12-31, replayed one week at a time with today’s PDB data.

- Released: 8,978; later added to mpstruc: 447; flagged at any likelihood: 444 (99.3%)

| Show if likelihood ≥ | Per week | Finds (recall) | Of those shown, in mpstruc |
|---|---|---|---|
| 0.00 | 137.1 | 99.3% | 12.3% |
| 0.10 | 100.2 | 99.3% | 16.9% |
| 0.20 | 76.8 | 99.1% | 21.9% |
| 0.30 | 61.3 | 96.6% | 26.8% |
| 0.40 | 50.8 | 91.9% | 30.8% |
| 0.41 ← default | 50.1 | 91.7% | 31.1% |
| 0.50 | 39.1 | 81.0% | 35.2% |
| 0.60 | 30.7 | 70.2% | 38.9% |
| 0.70 | 15.1 | 39.8% | 44.9% |
| 0.80 | 4.1 | 12.8% | 52.8% |
| 0.90 | 0.7 | 3.1% | 73.7% |

- Placement type right: 77.7% of all 444; 82.1% of the 420 the automation decides (5.4% go to curators with no default)
- Placement fully right (type + subfamily + parent/owner): 61.0% of all; 64.5% of those it decides
- Subfamily right first time: 92.8%; in top 3: 95.5% (of all 444; a subfamily was suggested for 99.5%)

| Rule | Structures | Type right | Fully right | Subfamily right | Parent/owner right |
|---|---|---|---|---|---|
| new-protein | 50 | 66.0% | 48.0% | 80.0% | — |
| no-evidence | 2 | 0.0% | 0.0% | None% | — |
| partial-protein | 22 | 0.0% | 0.0% | 86.4% | 58.8% |
| same-paper | 23 | 78.3% | 60.9% | 100.0% | 70.0% |
| same-paper-batch | 257 | 87.9% | 65.8% | 93.4% | 74.7% |
| same-protein | 90 | 75.6% | 71.1% | 100.0% | 73.6% |

Confusion (true → predicted): daughter→daughter: 68, daughter→parent: 4, daughter→related: 21, daughter→review: 10, parent→daughter: 3, parent→parent: 33, parent→related: 15, parent→review: 7, related→daughter: 19, related→parent: 13, related→related: 244, related→review: 7

## Tuning period, release week (in-sample: the model was fitted on these)

Structures released 2025-01-01 to 2025-06-30, replayed one week at a time as each looked in its release week.

- Released: 8,622; later added to mpstruc: 367; flagged at any likelihood: 366 (99.7%)

| Show if likelihood ≥ | Per week | Finds (recall) | Of those shown, in mpstruc |
|---|---|---|---|
| 0.00 | 105.8 | 99.7% | 13.4% |
| 0.10 | 78.0 | 99.7% | 18.2% |
| 0.20 | 65.9 | 99.5% | 21.4% |
| 0.30 | 55.3 | 99.5% | 25.5% |
| 0.40 | 43.9 | 93.7% | 30.3% |
| 0.41 ← default | 43.3 | 93.7% | 30.7% |
| 0.50 | 33.3 | 83.4% | 35.5% |
| 0.60 | 25.8 | 74.4% | 40.9% |
| 0.70 | 14.3 | 52.6% | 52.3% |
| 0.80 | 6.0 | 30.8% | 73.4% |
| 0.90 | 2.0 | 12.8% | 88.7% |

- Placement type right: 81.7% of all 366; 84.0% of the 356 the automation decides (2.7% go to curators with no default)
- Placement fully right (type + subfamily + parent/owner): 67.2% of all; 69.1% of those it decides
- Subfamily right first time: 90.2%; in top 3: 92.6% (of all 366; a subfamily was suggested for 96.4%)

| Rule | Structures | Type right | Fully right | Subfamily right | Parent/owner right |
|---|---|---|---|---|---|
| new-protein | 39 | 66.7% | 56.4% | 79.5% | — |
| no-evidence | 3 | 0.0% | 0.0% | None% | — |
| partial-protein | 7 | 0.0% | 0.0% | 85.7% | 80.0% |
| same-paper | 19 | 73.7% | 47.4% | 100.0% | 62.5% |
| same-paper-batch | 226 | 92.0% | 73.0% | 93.5% | 80.9% |
| same-protein | 72 | 70.8% | 69.4% | 100.0% | 74.6% |

Confusion (true → predicted): daughter→daughter: 51, daughter→parent: 2, daughter→related: 14, daughter→review: 5, parent→daughter: 5, parent→parent: 26, parent→related: 9, parent→review: 5, related→daughter: 16, related→parent: 11, related→related: 222

## Likelihood model (fitted on the tuning period)

| Signal | Weight |
|---|---|
| UniProt: transmembrane protein | +1.10 |
| Same protein as an mpstruc entry | +0.65 |
| Cryo-EM | +0.49 |
| Same paper as another candidate this week | +0.48 |
| Same paper as an mpstruc entry | +0.47 |
| UniProt: membrane protein | +0.45 |
| NMR | +0.02 |
| X-ray | -0.83 |

Title words that raise the likelihood most: inward-facing (+1.43), agonist (+1.09), coli (+0.93), lmng (+0.88), lipid (+0.85), methyl (+0.76), human (+0.75), scfv16 (+0.73), class (+0.71), msba (+0.70), intermediate (+0.62), obtained (+0.61), xpr1 (+0.58), inhibitor (+0.57), protein (+0.49)

Title words that lower it most: sars-cov-2 (-1.72), spike (-1.23), virus (-1.00), compound (-0.99), assembly (-0.91), dimer (-0.90), ca2 (-0.87), map (-0.83), nanodisc (-0.82), rbd (-0.78), protease (-0.71), glycoprotein (-0.70), supercomplex (-0.68), variant (-0.67), trimer (-0.67)

## Confidence shown for each rule (fully right, tuning period)

- parent: 56%
- partial: 0%
- related: 47%
- related-batch: 73%
- daughter: 69%
