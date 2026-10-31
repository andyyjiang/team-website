# Automation backtest

The automation is run as if it were the start of each period, knowing only the mpstruc entries released before then, and its proposals are compared with what the curators actually did.

## Tuning period: structures released 2025-01-01 to 2025-06-30

- Released in the PDB: 8,622
- Of these, later added to mpstruc: 367
- **Found by the automation: 328 (89.4%)**
- Proposed but not (yet) in mpstruc: 1,545 (not necessarily wrong: mpstruc may not have added them yet)
- **Placement type right: 65.5%** (of 328)
- Subfamily right (first choice): 86.7%; in top 3: 96.5%

| Rule | Structures | Type right | Subfamily right | Parent/owner right |
|---|---|---|---|---|
| new-protein | 23 | 43.5% | 78.3% | — |
| no-evidence | 2 | 0.0% | None% | 0.0% |
| partial-protein | 19 | 0.0% | 63.2% | 21.4% |
| same-paper | 5 | 80.0% | 100.0% | 80.0% |
| same-paper-batch | 229 | 76.0% | 86.3% | 21.1% |
| same-protein | 50 | 54.0% | 100.0% | 52.0% |

Confusion (true → predicted): daughter→daughter: 27, daughter→parent: 3, daughter→related: 36, daughter→review: 6, parent→parent: 10, parent→related: 20, parent→review: 5, related→daughter: 23, related→parent: 10, related→related: 178, related→review: 10

## Test period (untouched by tuning): structures released 2025-07-01 to 2025-12-31

- Released in the PDB: 8,978
- Of these, later added to mpstruc: 447
- **Found by the automation: 356 (79.6%)**
- Proposed but not (yet) in mpstruc: 1,391 (not necessarily wrong: mpstruc may not have added them yet)
- **Placement type right: 58.4%** (of 356)
- Subfamily right (first choice): 89.6%; in top 3: 98.0%

| Rule | Structures | Type right | Subfamily right | Parent/owner right |
|---|---|---|---|---|
| new-protein | 12 | 50.0% | 66.7% | — |
| partial-protein | 41 | 0.0% | 63.4% | 27.0% |
| same-paper | 7 | 85.7% | 100.0% | 85.7% |
| same-paper-batch | 236 | 72.0% | 92.4% | 19.7% |
| same-protein | 60 | 43.3% | 100.0% | 43.1% |

Confusion (true → predicted): daughter→daughter: 26, daughter→related: 54, daughter→review: 18, parent→daughter: 2, parent→parent: 6, parent→related: 13, parent→review: 4, related→daughter: 32, related→parent: 6, related→related: 176, related→review: 19

## Confidence used for each rule (from the tuning period)

- parent: 34%
- partial: 0%
- related: 80%
- related-batch: 66%
- daughter: 54%
