# Counting barrel strands: expected vs counted

Strands of the membrane beta-barrel counted from C-alpha coordinates (`mpstruc/barrel.py`), against the
number in the literature, for official mpstruc entries. Regenerate with `python -m mpstruc.barrel validate`.

**Exactly right: 56 of 57 (98%)**; within one strand: 56 of 57. In the two subfamilies with a curator rule: 27 of 27.

"Counted" is geometry only, over every chain. "Read as a membrane barrel" is what the membrane score and the
placement rules do with it (`membrane.counted_barrels`): only chains of proteins UniProt places in the outer
membrane are looked at, and the ring must lie in the protein's membrane part.

| PDB | Protein | In mpstruc under | Expected | Counted | Barrel | | Read as a membrane barrel |
|---|---|---|---|---|---|---|---|
| 1UYN | NalP translocator domain | G02-S006 | 12 | 12 | closed | ok | yes |
| 2QOM | EspP | G02-S006 | 12 | 12 | closed | ok | yes |
| 3SLJ | EspP, pre-cleavage | G02-S006 | 12 | 12 | closed | ok | yes |
| 3KVN | EstA | G02-S006 | 12 | 12 | closed | ok | yes |
| 3AEH | Hbp | G02-S006 | 12 | 12 | closed | ok | yes |
| 3QQ2 | BrkA | G02-S006 | 12 | 12 | closed | ok | yes |
| 4MEE | AIDA-I | G02-S006 | 12 | 12 | closed | ok | yes |
| 4E1S | Intimin | G02-S006 | 12 | 12 | closed | ok | yes |
| 4E1T | Invasin | G02-S006 | 12 | 12 | closed | ok | only when suggested under a ruled subfamily (UniProt names no outer-membrane protein) |
| 2GR8 | Hia (trimeric, 4 per chain) | G02-S006 | 12 | 12 | closed, 3 chains | ok | yes |
| 2LME | YadA (trimeric, 4 per chain) | G02-S006 | 12 | 12 | closed, 3 chains | ok | yes |
| 4C00 | TamA | G02-S006 | 16 | 16 | closed | ok | yes |
| 6WIL | CdiB (a TpsB transporter) | G02-S006 | 16 | 16 | closed | ok | yes |
| 6WIM | CdiB (a TpsB transporter) | G02-S006 | 16 | 16 | closed | ok | yes |
| 4QKY | FhaC | G02-S007 | 16 | 16 | closed | ok | yes |
| 3NJT | FhaC R450A | G02-S007 | 16 | 16 | closed | ok | yes |
| 4K3B | BamA | G02-S007 | 16 | 16 | closed | ok | yes |
| 4K3C | BamA | G02-S007 | 16 | 16 | closed | ok | yes |
| 4N75 | BamA barrel | G02-S007 | 16 | 16 | closed | ok | yes |
| 4C4V | BamA | G02-S007 | 16 | 16 | closed | ok | yes |
| 6QGW | BamA with nanobody | G02-S007 | 16 | 16 | closed | ok | yes |
| 5EKQ | BamACDE | G02-S007 | 16 | 16 | open seam | ok | yes |
| 5D0O | BamABCDE | G02-S007 | 16 | 16 | closed | ok | yes |
| 5LJO | BamABCDE (related code) | G02-S007 | 16 | 16 | open seam | ok | yes |
| 6V05 | BamABCDE with substrate | G02-S007 | 16 | 16 | open seam | ok | yes |
| 8BWC | BamABCDE | G02-S007 | 16 | 16 | open seam | ok | yes |
| 8US1 | TamA | G02-S007 | 16 | 16 | open seam | ok | yes |
| 1BXW | OmpA | control | 8 | 8 | closed | ok | yes |
| 1QJP | OmpA | control | 8 | 8 | closed | ok | yes |
| 1THQ | PagP | control | 8 | 8 | closed | ok | yes |
| 1QJ8 | OmpX | control | 8 | 8 | closed | ok | yes |
| 1K24 | OpcA | control | 10 | 10 | closed | ok | no |
| 1I78 | OmpT | control | 10 | 10 | closed | ok | yes |
| 1QD6 | OmpLA | control | 12 | 12 | closed | ok | yes |
| 1TLY | Tsx | control | 12 | 12 | closed | ok | yes |
| 3DWO | FadL | control | 14 | 14 | closed | ok | yes |
| 1T16 | FadL | control | 14 | 14 | closed | ok | yes |
| 2OMF | OmpF | control | 16 | 16 | closed | ok | yes |
| 1OPF | OmpF | control | 16 | 16 | closed | ok | yes |
| 2POR | porin | control | 16 | 16 | closed | ok | yes |
| 1PRN | porin | control | 16 | 16 | closed | ok | yes |
| 1A0S | ScrY sucrose porin | control | 18 | 18 | closed | ok | yes |
| 1MAL | LamB maltoporin | control | 18 | 18 | closed | ok | yes |
| 2O4V | OprP | control | 16 | 16 | closed | ok | yes |
| 3SY7 | OccD1 (OprD) | control | 18 | 18 | closed | ok | yes |
| 1FEP | FepA (TonB-dependent) | control | 22 | 22 | closed | ok | yes |
| 2FCP | FhuA (TonB-dependent) | control | 22 | 22 | closed | ok | yes |
| 1KMO | FecA (TonB-dependent) | control | 22 | 22 | closed | ok | yes |
| 1NQE | BtuB (TonB-dependent) | control | 22 | 22 | closed | ok | yes |
| 4Q35 | LptD | control | 26 | 26 | closed | ok | yes |
| 5IV9 | LptD | control | 26 | 26 | closed | ok | yes |
| 2J1N | OmpC | control | 16 | 16 | closed | ok | yes |
| 3EMN | VDAC1 (mitochondrial) | control | 19 | 19 | closed | ok | yes |
| 2JK4 | VDAC1 (mitochondrial) | control | 19 | 19 | closed | ok | yes |
| 1EK9 | TolC (trimer, 4 per chain) | control | 12 | 12 | closed, 3 chains | ok | yes |
| 7AHL | alpha-hemolysin (heptamer, 2 per chain) | control | 14 | 36 | closed, 7 chains | **differs** | no |
| 1UUN | MspA (octamer, 2 per chain; two chains deposited) | control | 16 | 16 | closed, 8 chains | ok | yes |

## Barrels with an open seam

A barrel caught with a substrate in it or while it folds has no closed ring; its count is a lower bound and is
shown as "at least N". Exactly right: **14 of 16** open barrels (those in the table above
and the Omp85 barrels below, 16 strands each).

| PDB | Protein | In mpstruc under | Expected | Counted | Barrel | | Read as a membrane barrel |
|---|---|---|---|---|---|---|---|
| 7RI4 | BamA with a substrate (EspP) | open seam | 16 | 15 | open seam | **differs** | yes |
| 7TT3 | BamA with a substrate | open seam | 16 | 16 | open seam | ok | yes |
| 7TSZ | BamA with a substrate | open seam | 16 | 16 | open seam | ok | yes |
| 7RJ5 | BamA with a substrate | open seam | 16 | 16 | open seam | ok | yes |
| 9SOW | BamA with a substrate (new structure) | open seam | 16 | 16 | open seam | ok | yes |
| 6WUH | Sam50 with a substrate (Tom40) | open seam | 16 | 16 | open seam | ok | yes |
| 6WUJ | Sam50 with a substrate | open seam | 16 | 15 | open seam | **differs** | yes |
| 6WUT | Sam50 with a substrate | open seam | 16 | 16 | open seam | ok | yes |
| 6WUM | Sam50 with a substrate | open seam | 16 | 16 | open seam | ok | yes |
| 25WF | SAM50 (new structure) | open seam | 16 | 16 | open seam | ok | yes |
| 9TGQ | TamA with TamB (new structure) | open seam | 16 | 16 | open seam | ok | yes |

## Barrels and sheets that are not membrane barrels

The count is geometry only: soluble barrels close too. They are never read as membrane barrels because
the membrane score and the placement rules also require UniProt to place the protein in the outer membrane
(see `membrane.OM_FAMILIES`, `membrane.counted_barrels`).

| PDB | Protein | Counted (geometry) | Read as a membrane barrel |
|---|---|---|---|
| 1EMA | GFP | no barrel | no |
| 1STP | streptavidin | 8 (closed) | no |
| 1BBP | bilin-binding protein (lipocalin) | 8 (closed) | no |
| 8TIM | triosephosphate isomerase (TIM barrel) | 8 (closed) | no |
| 3PRW | BamB (beta-propeller) | no barrel | no |
| 3VE2 | TbpB (lipoprotein on the outer membrane) | 8 (closed) | no |
| 3KVD | fHbp (lipoprotein on the outer membrane) | 8 (closed) | no |

## Guard set: soluble fragments of outer-membrane proteins

UniProt's keywords describe the full-length protein, so a head, stalk, passenger, POTRA or periplasmic
fragment carries them too. A ring of strands in such a fragment is not read as a membrane barrel: it must lie
in the membrane part UniProt names (transmembrane strands, a barrel domain, the membrane anchor or
translocator domain), and a ring made by several chains that cannot be checked that way (UniProt names no
membrane part) needs 12 or more strands, a whole number per chain (`membrane.membrane_barrel`). A ring of fewer
than 8 strands is never a membrane barrel.

**Read as a membrane barrel: 0 of 32** (must be 0).

| PDB | Protein | Counted (geometry) | Read as a membrane barrel |
|---|---|---|---|
| 1S7M | Hia binding domain HiaBD1 (trimeric autotransporter head) | 8 (closed, 3 chains) | no |
| 3EMI | Hia 307-422 (trimeric autotransporter) | 6 (closed, 3 chains) | no |
| 3LA9 | BpaA head domain (trimeric autotransporter) | 8 (closed, 3 chains) | no |
| 4USX | BPSL2063 YadA-like head (trimeric autotransporter) | 9 (closed, 3 chains) | no |
| 3D9X | BadA head (trimeric autotransporter) | 15 (closed, 3 chains) | no |
| 3WPR | AtaA stalk fragment fused to GCN4 (trimeric autotransporter) | 17 (closed, 3 chains) | no |
| 2YO2 | SadA 255-358 fused to GCN4 (trimeric autotransporter) | 9 (closed, 3 chains) | no |
| 2YNZ | SadA 823-947 fused to GCN4 (trimeric autotransporter) | 7 (closed, 3 chains) | no |
| 2XQH | EibD immunoglobulin-binding fragment (trimeric autotransporter) | 6 (closed, 3 chains) | no |
| 2QCZ | BamA POTRA domains | no barrel | no |
| 2QDF | BamA POTRA domains | no barrel | no |
| 3EFC | BamA periplasmic domain | no barrel | no |
| 3OG5 | BamA POTRA 4-5 | no barrel | no |
| 4XGA | BamB with BamA POTRA 3-5 | no barrel | no |
| 3MC8 | Omp85 POTRA 1-3 (Anabaena) | no barrel | no |
| 2X8X | Omp85 N-terminal domain (Thermosynechococcus) | no barrel | no |
| 5UAY | Toc75 POTRA domains | no barrel | no |
| 1WXR | Hbp passenger domain (autotransporter) | 4 (closed) | no |
| 1DAB | pertactin P.69 passenger domain (autotransporter) | no barrel | no |
| 1F02 | intimin extracellular fragment | no barrel | no |
| 3BWU | FimD usher N-terminal domain with FimC-FimH | 4 (closed) | no |
| 1ZE3 | FimD usher N-terminal domain with FimC-FimH pilin domain | no barrel | no |
| 3L48 | PapC usher C-terminal domain | 4 (closed) | no |
| 1ZZV | FecA periplasmic signalling domain (TonB-dependent receptor) | no barrel | no |
| 2A02 | PupA periplasmic signalling domain (TonB-dependent receptor) | no barrel | no |
