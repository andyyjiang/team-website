# Automation backtest

Each week of a past period is replayed: the automation knows only the mpstruc entries released before that week and is scored against what mpstruc's curators eventually did. The likelihood model, the default cut-off and the per-rule confidences are fitted on Jan–Jun 2025 only (cut-off by cross-validation grouped by paper). The placement rules themselves were designed from the full 2 Oct 2026 snapshot (reports/placement_analysis.md), which includes the test period (about 4% of it).

## Summary: test period (Jul–Dec 2025, never used for fitting), at the default cut-off

| Data the automation has | Finds | Per week | Shown that mpstruc added | Subfamily 1st (top 3) | Type right* | Fully right* | No default |
|---|---|---|---|---|---|---|---|
| release week | 92.8% | 52.9 | 29.8% | 99.1% (99.4%) | 58.9% | 47.3% | 27.8% |
| once papers are linked | 98.0% | 69.4 | 24.0% | 93.2% (95.9%) | 82.1% | 64.5% | 5.4% |

\* of the structures the automation places itself. *Release week* is what curators see first. The weekly job re-checks the last 26 weeks every Thursday, so as papers are linked to a structure its suggestion moves toward the *once papers are linked* row. About half of new releases have no paper linked in their first weeks.

## Test period, release week

Structures released 2025-07-01 to 2025-12-31, replayed one week at a time as each looked in its release week.

- Released: 8,978; later added to mpstruc: 447; flagged at any likelihood: 442 (98.9%)

| Show if likelihood ≥ | Per week | Finds (recall) | Of those shown, in mpstruc |
|---|---|---|---|
| 0.00 | 106.1 | 98.9% | 15.9% |
| 0.10 | 80.1 | 98.9% | 21.0% |
| 0.20 | 64.0 | 98.4% | 26.2% |
| 0.30 | 55.9 | 94.4% | 28.7% |
| 0.35 ← default | 52.9 | 92.8% | 29.8% |
| 0.40 | 44.3 | 87.0% | 33.4% |
| 0.50 | 37.1 | 81.0% | 37.1% |
| 0.60 | 22.8 | 60.4% | 45.1% |
| 0.70 | 14.8 | 49.0% | 56.2% |
| 0.80 | 7.8 | 34.0% | 73.8% |
| 0.90 | 2.2 | 10.3% | 78.0% |

- Placement type right: 42.5% of all 442; 58.9% of the 319 the automation decides (27.8% go to curators with no default)
- Placement fully right (type + subfamily + parent/owner): 34.2% of all; 47.3% of those it decides
- Subfamily right first time: 99.1%; in top 3: 99.4%

| Rule | Structures | Type right | Fully right | Subfamily right | Parent/owner right |
|---|---|---|---|---|---|
| no-evidence | 78 | 0.0% | 0.0% | None% | 0.0% |
| partial-protein | 45 | 0.0% | 0.0% | 93.3% | 30.8% |
| same-paper | 7 | 85.7% | 57.1% | 100.0% | 57.1% |
| same-paper-batch | 120 | 87.5% | 61.7% | 100.0% | 79.8% |
| same-protein | 192 | 40.1% | 38.0% | 100.0% | 39.6% |

Confusion (true → predicted): daughter→daughter: 77, daughter→related: 10, daughter→review: 15, parent→daughter: 5, parent→related: 6, parent→review: 46, related→daughter: 110, related→related: 111, related→review: 62

## Test period, once papers are linked

Structures released 2025-07-01 to 2025-12-31, replayed one week at a time with today’s PDB data.

- Released: 8,978; later added to mpstruc: 447; flagged at any likelihood: 444 (99.3%)

| Show if likelihood ≥ | Per week | Finds (recall) | Of those shown, in mpstruc |
|---|---|---|---|
| 0.00 | 137.1 | 99.3% | 12.3% |
| 0.10 | 108.3 | 99.3% | 15.6% |
| 0.20 | 85.0 | 99.1% | 19.8% |
| 0.30 | 74.0 | 98.4% | 22.6% |
| 0.35 ← default | 69.4 | 98.0% | 24.0% |
| 0.40 | 60.7 | 96.6% | 27.1% |
| 0.50 | 53.5 | 93.3% | 29.7% |
| 0.60 | 43.3 | 88.4% | 34.7% |
| 0.70 | 32.5 | 76.5% | 40.0% |
| 0.80 | 23.2 | 62.2% | 45.6% |
| 0.90 | 6.6 | 20.8% | 53.8% |

- Placement type right: 77.7% of all 444; 82.1% of the 420 the automation decides (5.4% go to curators with no default)
- Placement fully right (type + subfamily + parent/owner): 61.0% of all; 64.5% of those it decides
- Subfamily right first time: 93.2%; in top 3: 95.9%

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
| 0.10 | 74.6 | 99.7% | 19.0% |
| 0.20 | 62.8 | 99.5% | 22.5% |
| 0.30 | 53.0 | 99.5% | 26.6% |
| 0.35 ← default | 49.9 | 99.2% | 28.2% |
| 0.40 | 41.7 | 94.0% | 32.0% |
| 0.50 | 34.7 | 88.3% | 36.1% |
| 0.60 | 20.6 | 70.6% | 48.7% |
| 0.70 | 13.5 | 54.0% | 56.6% |
| 0.80 | 8.7 | 40.9% | 66.4% |
| 0.90 | 3.5 | 19.6% | 80.0% |

- Placement type right: 45.6% of all 366; 60.7% of the 275 the automation decides (24.9% go to curators with no default)
- Placement fully right (type + subfamily + parent/owner): 31.1% of all; 41.5% of those it decides
- Subfamily right first time: 99.6%; in top 3: 100.0%

| Rule | Structures | Type right | Fully right | Subfamily right | Parent/owner right |
|---|---|---|---|---|---|
| no-evidence | 82 | 0.0% | 0.0% | None% | 0.0% |
| partial-protein | 9 | 0.0% | 0.0% | 88.9% | 71.4% |
| same-paper | 10 | 80.0% | 70.0% | 100.0% | 77.8% |
| same-paper-batch | 114 | 90.4% | 45.6% | 100.0% | 80.9% |
| same-protein | 151 | 37.1% | 36.4% | 100.0% | 38.9% |

Confusion (true → predicted): daughter→daughter: 56, daughter→related: 8, daughter→review: 8, parent→daughter: 7, parent→related: 5, parent→review: 33, related→daughter: 88, related→related: 111, related→review: 50

## Likelihood model (fitted on the tuning period)

| Signal | Weight |
|---|---|
| Same paper as another candidate this week | +1.26 |
| UniProt: transmembrane protein | +1.19 |
| Same paper as an mpstruc entry | +0.80 |
| Same protein as an mpstruc entry | +0.71 |
| Cryo-EM | +0.58 |
| UniProt: membrane protein | +0.54 |
| NMR | +0.05 |
| More than 3 protein chains | +0.03 |
| X-ray | -0.95 |

Title words that raise the likelihood most: inward-facing (+1.55), agonist (+1.20), methyl (+1.10), coli (+0.94), inhibitor (+0.84), msba (+0.83), human (+0.81), lipid (+0.77), lmng (+0.63), intermediate (+0.58), class (+0.56), scfv16 (+0.53), complexed (+0.50), mitochondrial (+0.49), protein (+0.49)

Title words that lower it most: sars-cov-2 (-1.70), spike (-1.12), virus (-1.09), compound (-0.97), ca2 (-0.91), assembly (-0.84), nanodisc (-0.80), rbd (-0.79), dimer (-0.77), supercomplex (-0.72), map (-0.70), fab (-0.68), antibody (-0.62), protease (-0.62), glycoprotein (-0.59)

## Confidence shown for each rule (fully right, tuning period)

- partial: 0%
- related: 70%
- related-batch: 46%
- daughter: 36%
