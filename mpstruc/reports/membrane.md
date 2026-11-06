# How likely is a structure to be a membrane protein?

A structure counts as a membrane protein if mpstruc lists it, if PDBTM or MemProtMD (which screened every PDB entry for transmembrane proteins) do, or if OPM classes it as transmembrane. Structures OPM classes only as monotopic/peripheral or peptide, with no other source, are left out. Signals include UniProt keywords and families (including the beta-barrel and outer-membrane keywords) and TCDB classes. The strongest signal is whether the part that crosses the membrane is resolved in the structure: UniProt transmembrane segments, or for beta-barrels a UniProt barrel domain (such as 'TBDR beta-barrel'). The list of barrel domain names was chosen from the training period and from published domain definitions, without looking at the test period.

Facts from the structure itself (what is resolved, detergents, lipids, the title) and from mpstruc are as they were on the release day. UniProt, TCDB and InterPro are read as they are today, and some of their annotations may have been added after (or because of) the structure. So the test-period numbers below may look somewhat better than they will be for brand-new entries.

It was trained on PDB releases 2023-07-01 to 2024-08-31 (calibrated on structures each fold had not seen) and tested on 2024-09-01 to 2025-03-31.

## Test period

- Structures: 9,583, of which membrane proteins: 1,187
- Ranking quality (AUC): 0.9853; Brier score: 0.0248
- Called membrane protein (≥50%): 1,324; of those, 83.7% are; they include 93.3% of all membrane proteins
- Of the 451 structures mpstruc added from Jul–Dec 2025 releases, 90.5% score at least 50% and 98.0% at least 10%

### Which kinds are found?

Share of each kind in the test period that scores at least 50%: every structure in each mpstruc group, and the membrane proteins UniProt places in a bacterial outer membrane.

| Kind | Structures | Scored 50% or more |
|---|---|---|
| Monotopic membrane proteins | 7 | 0.0% |
| Transmembrane proteins: alpha-helical | 474 | 93.0% |
| Transmembrane proteins: beta-barrel | 19 | 89.5% |
| outer membrane (UniProt) | 32 | 84.4% |

### Wrong calls: model errors or missing labels?

216 test structures labelled 'not a membrane protein' score 50% or more. 4 of them are proteins UniProt calls outer-membrane beta-barrels (cell outer membrane and transmembrane beta strand), which the label sources (mpstruc, PDBTM, MemProtMD, OPM) probably missed; the same protein is sometimes labelled a membrane protein in another structure. The other 212 are real mistakes. The labels were not changed, so the precision above counts both kinds as wrong. Outer-membrane barrels counted as wrong: 8RM3, 8RMI, 8S34, 9FVQ.

### Beta-barrels

Scores from this model for structures curators asked about: 38FI 90%; 38AQ 87%.

mpstruc beta-barrel structures in the test period that still score below 50% (2). Barrels whose UniProt entry names no barrel domain (only a plug, or nothing) are still often missed:

| PDB | Shown | Title |
|---|---|---|
| 9B8A | 1% | Cryo-EM structure of E. coli cellulose synthase subunit C |
| 9B8H | 1% | Cryo-EM structure of E. coli cellulose synthase subunit C with cellotetraose |

### Does the percentage mean what it says?

| Shown as | Structures | Average shown | Actually membrane proteins |
|---|---|---|---|
| 0-10% | 8,066 | 0.4% | 0.3% |
| 10-20% | 71 | 14.3% | 15.5% |
| 20-30% | 38 | 25.1% | 28.9% |
| 30-40% | 45 | 33.5% | 28.9% |
| 40-50% | 39 | 44.2% | 56.4% |
| 50-60% | 28 | 55.4% | 35.7% |
| 60-70% | 28 | 64.9% | 60.7% |
| 70-80% | 108 | 76.2% | 85.2% |
| 80-90% | 191 | 85.5% | 73.8% |
| 90-100% | 969 | 94.3% | 87.5% |

## Strongest signals

| Signal | Weight |
|---|---|
| tm_modelled | +2.49 |
| detergent | +1.33 |
| tm_missing | -1.33 |
| kw_tm | +1.22 |
| m:xray | -1.02 |
| tcdb_other | +0.82 |
| kw_tm_all | +0.74 |
| kw_KW-0963 | -0.67 |
| kw_membrane | +0.53 |
| family | +0.47 |
| many_chains | -0.38 |
| m:em | +0.34 |
| no_uniprot | -0.29 |
| kw_outer | +0.29 |
| tm_unmapped | +0.29 |
| lipid | +0.26 |
| kw_KW-0539 | +0.25 |
| tcdb_transporter | +0.25 |
| barrel_modelled | -0.16 |
| kw_KW-0964 | -0.12 |
| tcdb_barrel | +0.09 |
| kw_strand | +0.02 |
| same_protein | +0.01 |

Title words that raise it: receptor (+0.51), substrate (+0.40), channel (+0.34), full-length (+0.32), transporter (+0.32), region (+0.28), nanodiscs (+0.28), apo (+0.26), system (+0.26), lipid (+0.23), dimer (+0.21), outward-facing (+0.20), vibrio (+0.19), coli (+0.17), primer (+0.15)

Title words that lower it: variant (-0.43), fab (-0.33), composite (-0.30), domain (-0.30), beta (-0.28), intermediate (-0.26), protein (-0.26), sars-cov-2 (-0.25), fibril (-0.25), antibody (-0.25), iii (-0.24), dataset (-0.24), rat (-0.23), bacteriophage (-0.23), mouse (-0.23)

## Examples from the test period

| PDB | Shown | Membrane protein? | Why |
|---|---|---|---|
| 8JC9 | 96% | yes | 2 membrane-crossing segments are resolved in the structure; its protein family (IPR000066) is a membrane protein in 10 of 11 training structures; the sample contains detergent (D12) |
| 8YMN | 95% | yes | 11 membrane-crossing segments are resolved in the structure; its protein family (IPR003864) is a membrane protein in 21 of 21 training structures; UniProt lists a protein in it as crossing the membrane |
| 9EX4 | 95% | no | 4 membrane-crossing segments are resolved in the structure; its protein family (IPR036719) is a membrane protein in 94 of 103 training structures; UniProt lists a protein in it as crossing the membrane |
| 9ITY | 94% | no | 8 membrane-crossing segments are resolved in the structure; its protein family (IPR005953) is a membrane protein in 6 of 6 training structures; UniProt lists all its proteins as crossing the membrane |
| 8P15 | 90% | no | 7 membrane-crossing segments are resolved in the structure; its protein family (IPR017452) is a membrane protein in 363 of 388 training structures; UniProt lists a protein in it as crossing the membrane; but UniProt says a protein in it is in the cytoplasm |
| 8XNL | 78% | yes | 14 membrane-crossing segments are resolved in the structure; UniProt lists a protein in it as crossing the membrane; its protein family (IPR009346) is a membrane protein in 28 of 46 training structures |
| 8YJB | 23% | no | its protein family (IPR016024) is a membrane protein in 6 of 320 training structures; UniProt says a protein in it is in the cytoplasm, although 5 membrane-crossing segments are resolved in the structure |
| 9C5B | 4% | no | the part that crosses the membrane is not resolved, so only an outside part of the protein is in the structure; UniProt says a protein in it is in the cytoplasm, although UniProt lists a protein in it as crossing the membrane |
| 9BGM | 1% | no | UniProt does not list any of its proteins as crossing the membrane |
| 8J47 | 1% | no | the part that crosses the membrane is not resolved, so only an outside part of the protein is in the structure; its protein family (IPR011993) is a membrane protein in 6 of 247 training structures; UniProt says a protein in it is in the cytoplasm, although UniProt lists a protein in it as crossing the membrane |
| 9B04 | 1% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR014746) is a membrane protein in 0 of 20 training structures, although UniProt places a protein in it in a membrane |
| 8XKR | 0% | no | its protein family (IPR006212) is a membrane protein in 0 of 68 training structures; the part that crosses the membrane is not resolved, so only an outside part of the protein is in the structure; UniProt says a protein in it is secreted, although UniProt lists a protein in it as crossing the membrane |
| 8V5S | 0% | no | the part that crosses the membrane is not resolved, so only an outside part of the protein is in the structure; the title mentions 'fab', although UniProt lists a protein in it as crossing the membrane |
| 9CIW | 0% | no | UniProt does not list any of its proteins as crossing the membrane |
| 8XE8 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm |
| 9FZ8 | 0% | no | the part that crosses the membrane is not resolved, so only an outside part of the protein is in the structure; its protein family (IPR012338) is a membrane protein in 2 of 150 training structures, although UniProt lists a protein in it as crossing the membrane |
| 8S4K | 0% | no | none of its proteins is linked to UniProt yet |
| 9MRB | 0% | no | none of its proteins is linked to UniProt yet |
| 9DPF | 0% | no | its protein family (IPR009003) is a membrane protein in 9 of 1443 training structures; the part that crosses the membrane is not resolved, so only an outside part of the protein is in the structure; UniProt says a protein in it is secreted, although UniProt lists a protein in it as crossing the membrane |
| 9FVF | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR023631) is a membrane protein in 0 of 4 training structures |
| 8X32 | 0% | no | its protein family (IPR036388) is a membrane protein in 0 of 479 training structures; UniProt does not list any of its proteins as crossing the membrane, although UniProt says a protein in it is in the nucleus |
| 8UTI | 0% | no | its protein family (IPR036388) is a membrane protein in 0 of 479 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm, although TCDB lists a protein in it (TC 8.A.92.1.3) |
| 8KDZ | 0% | no | its protein family (IPR007094) is a membrane protein in 1 of 1429 training structures; the part that crosses the membrane is not resolved, so only an outside part of the protein is in the structure; UniProt says a protein in it is secreted, although UniProt lists a protein in it as crossing the membrane |
| 9H3T | 0% | no | its protein family (IPR009000) is a membrane protein in 0 of 340 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 9FAH | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR009013) is a membrane protein in 0 of 11 training structures |
| 9CBN | 0% | no | its protein family (IPR029053) is a membrane protein in 0 of 256 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 8ZQ2 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR002073) is a membrane protein in 0 of 16 training structures; UniProt says a protein in it is in the cytoplasm, although UniProt places a protein in it in a membrane |
| 8RSZ | 0% | no | its protein family (IPR001926) is a membrane protein in 0 of 32 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 9NB7 | 0% | no | its protein family (IPR000169) is a membrane protein in 0 of 32 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 9MM2 | 0% | no | its protein family (IPR012334) is a membrane protein in 0 of 26 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is secreted |
| 7GA8 | 0% | no | its protein family (IPR029071) is a membrane protein in 0 of 396 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm, although UniProt places a protein in it in a membrane |
| 7HAJ | 0% | no | its protein family (IPR020568) is a membrane protein in 0 of 385 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm, although UniProt places a protein in it in a membrane |
| 7HNV | 0% | no | its protein family (IPR001841) is a membrane protein in 0 of 192 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm, although UniProt says a protein in it is in the nucleus |
| 8JZQ | 0% | no | its protein family (IPR002213) is a membrane protein in 0 of 44 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 8S39 | 0% | no | its protein family (IPR046346) is a membrane protein in 0 of 37 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 8WTL | 0% | no | its protein family (IPR011545) is a membrane protein in 0 of 126 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 8YYJ | 0% | no | its protein family (IPR036866) is a membrane protein in 0 of 121 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm |
| 9C5D | 0% | no | its protein family (IPR002130) is a membrane protein in 0 of 40 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm |
| 9F15 | 0% | no | its protein family (IPR001148) is a membrane protein in 0 of 78 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm, although UniProt places a protein in it in a membrane |
| 9GVY | 0% | no | its protein family (IPR000566) is a membrane protein in 0 of 31 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm, although UniProt places a protein in it in a membrane |
| 9MTO | 0% | no | its protein family (IPR016181) is a membrane protein in 0 of 36 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm |
