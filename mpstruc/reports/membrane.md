# How likely is a structure to be a membrane protein?

A structure counts as a membrane protein if mpstruc lists it, if PDBTM or MemProtMD (which screened every PDB entry for transmembrane proteins) do, or if OPM classes it as transmembrane. Structures OPM classes only as monotopic/peripheral or peptide, with no other source, are left out. Signals include UniProt keywords and families and TCDB transporter classes. The model uses only facts known on the day a structure is released. It was trained on PDB releases 2023-07-01 to 2024-08-31 (calibrated on structures each fold had not seen) and tested on 2024-09-01 to 2025-03-31.

## Test period

- Structures: 9,583, of which membrane proteins: 1,187
- Ranking quality (AUC): 0.9832; Brier score: 0.0322
- Called membrane protein (≥50%): 1,317; of those, 79.8% are; they include 88.5% of all membrane proteins
- Of the 447 structures mpstruc added from Jul–Dec 2025 releases, 91.7% score at least 50% and 97.1% at least 10%

### Does the percentage mean what it says?

| Shown as | Structures | Average shown | Actually membrane proteins |
|---|---|---|---|
| 0-10% | 7,969 | 0.5% | 0.5% |
| 10-20% | 156 | 13.9% | 23.7% |
| 20-30% | 55 | 25.1% | 23.6% |
| 30-40% | 37 | 34.2% | 18.9% |
| 40-50% | 49 | 44.5% | 73.5% |
| 50-60% | 35 | 55.0% | 48.6% |
| 60-70% | 85 | 65.1% | 67.1% |
| 70-80% | 132 | 76.2% | 60.6% |
| 80-90% | 324 | 86.6% | 79.9% |
| 90-100% | 741 | 92.4% | 86.1% |

## Strongest signals

| Signal | Weight |
|---|---|
| kw_tm | +2.06 |
| detergent | +1.71 |
| m:xray | -1.36 |
| tcdb_other | +0.78 |
| kw_membrane | +0.67 |
| kw_tm_all | +0.62 |
| lipid | +0.62 |
| m:em | +0.62 |
| kw_KW-0963 | -0.57 |
| family | +0.57 |
| tcdb_transporter | +0.55 |
| kw_KW-0964 | -0.47 |
| no_uniprot | -0.30 |
| kw_KW-0539 | +0.25 |
| same_protein | +0.22 |
| m:nmr | +0.04 |
| many_chains | -0.04 |

Title words that raise it: transporter (+0.55), channel (+0.54), receptor (+0.52), substrate (+0.44), system (+0.31), coli (+0.31), nanodiscs (+0.30), apo (+0.28), outward-facing (+0.26), vibrio (+0.25), tuberculosis (+0.25), full-length (+0.24), membrane (+0.24), mycobacterium (+0.24), lipid (+0.23)

Title words that lower it: domain (-0.90), fab (-0.81), sars-cov-2 (-0.69), rat (-0.50), antibody (-0.48), isoform (-0.47), spike (-0.44), flip (-0.44), beta (-0.40), ampa (-0.39), local (-0.38), glua2 (-0.37), tarp (-0.36), protein (-0.34), composite (-0.33)

## Examples from the test period

| PDB | Shown | Membrane protein? | Why |
|---|---|---|---|
| 8KDF | 93% | yes | UniProt lists all its proteins as crossing the membrane; its protein family (IPR002293) is a membrane protein in 8 of 8 training structures; lipids are bound (CLR, LBN) |
| 9JUL | 93% | yes | its protein family (IPR011527) is a membrane protein in 93 of 93 training structures; UniProt lists a protein in it as crossing the membrane; TCDB classifies a protein in it as a channel or transporter (TC 3.A.1.201.5) |
| 9CTF | 93% | yes | its protein family (IPR011527) is a membrane protein in 93 of 93 training structures; UniProt lists a protein in it as crossing the membrane; TCDB classifies a protein in it as a channel or transporter (TC 3.A.1.201.1); but UniProt says a protein in it is in the cytoplasm |
| 8ZPM | 90% | no | its protein family (IPR017452) is a membrane protein in 363 of 388 training structures; UniProt lists a protein in it as crossing the membrane; TCDB lists a protein in it (TC 9.A.14.13.17) |
| 8WPN | 86% | yes | its protein family (IPR005821) is a membrane protein in 194 of 216 training structures; UniProt lists a protein in it as crossing the membrane; the sample contains detergent (Y01) |
| 9DYN | 70% | yes | UniProt lists a protein in it as crossing the membrane; TCDB classifies a protein in it as a channel or transporter (TC 1.A.46.1.2); mpstruc already lists this protein (e.g. 8D1E) |
| 9BVP | 26% | yes | its protein family (IPR011051) is a membrane protein in 0 of 22 training structures; but UniProt lists all its proteins as crossing the membrane |
| 9FX3 | 8% | no | UniProt does not list any of its proteins as crossing the membrane; but UniProt places a protein in it in a membrane |
| 9BGZ | 5% | no | none of its proteins is linked to UniProt yet; but the title mentions 'transporter' |
| 9HBW | 1% | no | UniProt does not list any of its proteins as crossing the membrane |
| 8QHZ | 1% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR007157) is a membrane protein in 0 of 11 training structures; but UniProt places a protein in it in a membrane |
| 9KBQ | 1% | no | UniProt does not list any of its proteins as crossing the membrane; but its protein family (IPR002196) is a membrane protein in 7 of 9 training structures |
| 9C3A | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR005564) is a membrane protein in 0 of 4 training structures |
| 8UDY | 0% | no | its protein family (IPR001205) is a membrane protein in 0 of 1186 training structures; the title mentions 'sars-cov-2'; but UniProt lists a protein in it as crossing the membrane |
| 9AR1 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; but its protein family (IPR001509) is a membrane protein in 17 of 64 training structures |
| 9CFU | 0% | no | its protein family (IPR043129) is a membrane protein in 0 of 123 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm; but UniProt places a protein in it in a membrane |
| 8QGA | 0% | no | UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm; but its protein family (IPR016064) is a membrane protein in 1 of 2 training structures |
| 8Z1G | 0% | no | its protein family (IPR036866) is a membrane protein in 0 of 121 training structures; UniProt does not list any of its proteins as crossing the membrane; but UniProt says a protein in it is in the nucleus |
| 9CQ9 | 0% | no | its protein family (IPR029063) is a membrane protein in 1 of 1382 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 9JJG | 0% | no | its protein family (IPR001205) is a membrane protein in 0 of 1186 training structures; UniProt does not list any of its proteins as crossing the membrane |
| 9E7Z | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR016184) is a membrane protein in 0 of 21 training structures; the title mentions 'fab' |
| 9CM2 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR009003) is a membrane protein in 9 of 1443 training structures; UniProt says a protein in it is secreted |
| 9JFF | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR004838) is a membrane protein in 0 of 3 training structures |
| 8RN9 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR007099) is a membrane protein in 0 of 42 training structures; the title mentions 'local' |
| 9MGX | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR000845) is a membrane protein in 0 of 7 training structures |
| 8YD1 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR029045) is a membrane protein in 0 of 63 training structures; UniProt says a protein in it is in the cytoplasm |
| 9AYP | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR014719) is a membrane protein in 0 of 14 training structures |
| 9A4K | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR029045) is a membrane protein in 0 of 63 training structures; UniProt says a protein in it is in the cytoplasm |
| 8V57 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR000169) is a membrane protein in 0 of 32 training structures; UniProt says a protein in it is secreted; but UniProt places a protein in it in a membrane |
| 7GUL | 0% | no | its protein family (IPR013087) is a membrane protein in 0 of 114 training structures; UniProt does not list any of its proteins as crossing the membrane; the title mentions 'domain'; but UniProt says a protein in it is in the nucleus |
| 7HBB | 0% | no | its protein family (IPR020568) is a membrane protein in 0 of 385 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm; but UniProt places a protein in it in a membrane |
| 7HUJ | 0% | no | its protein family (IPR001205) is a membrane protein in 0 of 1186 training structures; UniProt does not list any of its proteins as crossing the membrane; but UniProt places a protein in it in a membrane |
| 8P3H | 0% | no | its protein family (IPR035892) is a membrane protein in 0 of 107 training structures; UniProt does not list any of its proteins as crossing the membrane; UniProt says a protein in it is in the cytoplasm; but TCDB lists a protein in it (TC 9.A.57.1.4) |
| 8RW8 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR013767) is a membrane protein in 0 of 59 training structures; the title mentions 'domain'; but UniProt says a protein in it is in the nucleus |
| 8W13 | 0% | no | its protein family (IPR036388) is a membrane protein in 0 of 479 training structures; UniProt does not list any of its proteins as crossing the membrane; the title mentions 'domain'; but UniProt says a protein in it is in the nucleus |
| 8Y23 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR005468) is a membrane protein in 0 of 17 training structures; UniProt says a protein in it is secreted |
| 9AVO | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR019931) is a membrane protein in 0 of 46 training structures; the title mentions 'fab'; but UniProt says a protein in it is in the nucleus |
| 9D7L | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR016161) is a membrane protein in 0 of 37 training structures; the title mentions 'domain' |
| 9FD9 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR013785) is a membrane protein in 1 of 142 training structures; UniProt says a protein in it is in the cytoplasm |
| 9H7I | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR000895) is a membrane protein in 0 of 47 training structures; UniProt says a protein in it is in the cytoplasm |
| 9MY5 | 0% | no | UniProt does not list any of its proteins as crossing the membrane; its protein family (IPR025110) is a membrane protein in 0 of 58 training structures; the title mentions 'domain' |
