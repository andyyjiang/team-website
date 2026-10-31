# How likely is a structure to be a membrane protein?

A structure counts as a membrane protein if mpstruc lists it, if PDBTM or MemProtMD (which screened every PDB entry for transmembrane proteins) do, or if OPM classes it as transmembrane. Structures OPM classes only as monotopic/peripheral or peptide, with no other source, are left out. Signals include UniProt keywords and families and TCDB transporter classes. The model uses only facts known on the day a structure is released. It was trained on PDB releases 2023-07-01 to 2024-08-31 (calibrated on structures each fold had not seen) and tested on 2024-09-01 to 2025-03-31.

## Test period

- Structures: 9,583, of which membrane proteins: 1,187
- Ranking quality (AUC): 0.9852; Brier score: 0.0254
- Called membrane protein (≥50%): 1,301; of those, 84.1% are; they include 92.2% of all membrane proteins
- Of the 447 structures mpstruc added from Jul–Dec 2025 releases, 91.5% score at least 50% and 98.0% at least 10%

### Does the percentage mean what it says?

| Shown as | Structures | Average shown | Actually membrane proteins |
|---|---|---|---|
| 0-10% | 8,026 | 0.4% | 0.3% |
| 10-20% | 106 | 13.7% | 11.3% |
| 20-30% | 51 | 24.7% | 25.5% |
| 30-40% | 62 | 35.0% | 41.9% |
| 40-50% | 37 | 44.8% | 51.4% |
| 50-60% | 27 | 55.4% | 48.1% |
| 60-70% | 39 | 67.2% | 74.4% |
| 70-80% | 61 | 76.2% | 78.7% |
| 80-90% | 182 | 85.6% | 75.8% |
| 90-100% | 992 | 94.2% | 87.3% |

## Strongest signals

| Signal | Weight |
|---|---|
| tm_modelled | +2.31 |
| detergent | +1.38 |
| kw_tm | +1.37 |
| tm_missing | -1.37 |
| m:xray | -1.03 |
| tcdb_other | +0.81 |
| kw_membrane | +0.61 |
| kw_KW-0963 | -0.60 |
| kw_tm_all | +0.50 |
| family | +0.44 |
| tcdb_transporter | +0.38 |
| m:em | +0.36 |
| many_chains | -0.34 |
| no_uniprot | -0.29 |
| lipid | +0.21 |
| kw_KW-0539 | +0.20 |
| kw_KW-0964 | -0.15 |
| same_protein | +0.10 |
| m:nmr | +0.02 |

Title words that raise it: receptor (+0.55), substrate (+0.40), transporter (+0.39), channel (+0.35), full-length (+0.32), apo (+0.28), system (+0.28), region (+0.28), nanodiscs (+0.24), vibrio (+0.23), lipid (+0.19), dimer (+0.19), outward-facing (+0.18), coli (+0.18), tuberculosis (+0.16)

Title words that lower it: fab (-0.58), variant (-0.42), domain (-0.41), antibody (-0.38), hiv-1 (-0.34), beta (-0.31), composite (-0.30), env (-0.29), sars-cov-2 (-0.29), fibril (-0.27), protein (-0.26), iii (-0.25), dataset (-0.24), bacteriophage (-0.23), trimer (-0.22)

## Examples from the test period

| PDB | Shown | Membrane protein? | Why |
|---|---|---|---|
| 8IRJ | 95% | yes | 7 of its membrane-crossing segments are resolved in the structure; its protein family (IPR017978) is a membrane protein in 26 of 27 training structures; UniProt lists a protein in it as crossing the membrane |
| 9CZK | 95% | yes | 10 of its membrane-crossing segments are resolved in the structure; its protein family (IPR005821) is a membrane protein in 194 of 216 training structures; UniProt lists all its proteins as crossing the membrane |
| 9J4X | 95% | yes | 11 of its membrane-crossing segments are resolved in the structure; UniProt lists a protein in it as crossing the membrane; its protein family (IPR004331) is a membrane protein in 6 of 8 training structures |
| 8W8X | 94% | yes | 7 of its membrane-crossing segments are resolved in the structure; its protein family (IPR000367) is a membrane protein in 166 of 170 training structures; UniProt lists a protein in it as crossing the membrane; but UniProt says a protein in it is secreted |
| 8W7R | 91% | yes | 7 of its membrane-crossing segments are resolved in the structure; the sample contains detergent (OLC); UniProt lists a protein in it as crossing the membrane |
| 9F5Z | 77% | no | 9 of its membrane-crossing segments are resolved in the structure; UniProt lists a protein in it as crossing the membrane; its protein family (IPR036909) is a membrane protein in 39 of 67 training structures |
| 9B43 | 22% | yes | but UniProt lists a protein in it as crossing the membrane |
| 9FTI | 4% | no | UniProt does not list any of its proteins as crossing the membrane; but the title mentions 'channel' |
| 9HBZ | 2% | no | UniProt does not list any of its proteins as crossing the membrane |
| 9B4V | 1% | no | UniProt does not list any of its proteins as crossing the membrane |
| 9FLX | 1% | no | none of its proteins is linked to UniProt yet |
| 9MML | 1% | no | the title mentions 'fab'; none of its proteins is linked to UniProt yet |
| 9MU6 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR023213) is a membrane protein in 0 of 26 training structures |
| 8QYC | 0% | no | its protein family (IPR014001) is a membrane protein in 1 of 304 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 9IIA | 0% | no | UniProt does not list any of its proteins as crossing the membrane |
| 8UL5 | 0% | no | UniProt does not list any of its proteins as crossing the membrane |
| 9BFF | 0% | no | its protein family (IPR025110) is a membrane protein in 0 of 58 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 8WPJ | 0% | no | none of its proteins is linked to UniProt yet |
| 9AUS | 0% | no | none of its proteins is linked to UniProt yet |
| 8S3L | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR005103) is a membrane protein in 0 of 2 training structures; UniProt says a protein in it is secreted |
| 8W3B | 0% | no | its protein family (IPR001245) is a membrane protein in 1 of 252 training structures; its membrane-crossing segments are not resolved, so only an outside part of the protein is in the structure; UniProt says a protein in it is secreted; but UniProt lists a protein in it as crossing the membrane |
| 9B78 | 0% | no | its protein family (IPR043129) is a membrane protein in 0 of 123 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm |
| 9H3N | 0% | no | its protein family (IPR009000) is a membrane protein in 0 of 340 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 8YO3 | 0% | no | its protein family (IPR020568) is a membrane protein in 0 of 385 training structures; UniProt does not list any of its proteins as crossing the membrane; the title mentions 'domain' |
| 8ZMK | 0% | no | its protein family (IPR014729) is a membrane protein in 0 of 89 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm |
| 8ZXC | 0% | no | its protein family (IPR001487) is a membrane protein in 0 of 90 training structures; UniProt does not list any of its proteins as crossing the membrane; the title mentions 'domain'; but UniProt says a protein in it is in the nucleus |
| 9CGA | 0% | no | its protein family (IPR029062) is a membrane protein in 0 of 41 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm; but UniProt places a protein in it in a membrane |
| 9GO5 | 0% | no | its protein family (IPR043129) is a membrane protein in 0 of 123 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm |
| 7H6H | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR009003) is a membrane protein in 9 of 1443 training structures; UniProt says a protein in it is in the cytoplasm; but UniProt places a protein in it in a membrane |
| 7G37 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR017850) is a membrane protein in 0 of 19 training structures; UniProt says a protein in it is secreted |
| 7GVN | 0% | no | its protein family (IPR013087) is a membrane protein in 0 of 114 training structures; UniProt does not list any of its proteins as crossing the membrane; the title mentions 'domain'; but UniProt says a protein in it is in the nucleus |
| 7HGC | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR000313) is a membrane protein in 0 of 19 training structures; UniProt says a protein in it is in the cytoplasm; but UniProt says a protein in it is in the nucleus |
| 7HUP | 0% | no | its protein family (IPR001205) is a membrane protein in 0 of 1186 training structures; UniProt does not list any of its proteins as crossing the membrane; its membrane-crossing segments are not resolved, so only an outside part of the protein is in the structure; but UniProt places a protein in it in a membrane |
| 8Q2F | 0% | no | its protein family (IPR001128) is a membrane protein in 0 of 95 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm |
| 8S7T | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR006094) is a membrane protein in 0 of 29 training structures |
| 8WWQ | 0% | no | its protein family (IPR029063) is a membrane protein in 1 of 1382 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 8Z29 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR046346) is a membrane protein in 0 of 37 training structures |
| 9C4P | 0% | no | its protein family (IPR029035) is a membrane protein in 0 of 80 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 9ETR | 0% | no | its protein family (IPR001357) is a membrane protein in 0 of 66 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm; but UniProt says a protein in it is in the nucleus |
| 9GU5 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR047575) is a membrane protein in 0 of 34 training structures; UniProt says a protein in it is in the cytoplasm |
| 9MTO | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR016181) is a membrane protein in 0 of 36 training structures; UniProt says a protein in it is in the cytoplasm |
