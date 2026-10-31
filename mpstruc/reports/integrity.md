# mpstruc data-integrity report

Snapshot `Fri Oct 02 16:25:59 PDT 2026`, SHA-256 `4174f67bb5685afe0c0cf046f0c5bc5ceaa53162f8fd9568cf679f5187ea9f7a`.

Findings are for curator review. Nothing in the data has been changed.

| Check | Severity | Count |
|---|---|---|
| [PDB ID points to an unrelated structure](#pdb_id_points_elsewhere) | error 7 | 7 |
| [PDB IDs the PDB does not recognise](#pdb_unknown) | error 3 | 3 |
| [Malformed DOIs](#bad_doi) | warning 2 | 2 |
| [Malformed entity spellings](#broken_entity) | warning 3 | 3 |
| [Primary citation differs from the PDB entry's citation](#citation_differs_from_pdb) | warning 22, info 2 | 24 |
| [Citation title matches the PDB's, but its PubMed ID/DOI do not](#citation_ids_wrong) | warning 2 | 2 |
| [mpstruc cites a preprint that has since been published](#citation_preprint) | warning 2 | 2 |
| [Description cites a structure missing from the related list](#described_not_related) | warning 14 | 14 |
| [PubMed ID matches the PDB but the DOI does not](#doi_conflict) | warning 3, info 5 | 8 |
| [Unexpected or malformed tags](#malformed_markup) | warning 2 | 2 |
| [Empty name, species or taxonomic domain](#missing_field) | warning 14 | 14 |
| [PDB IDs in the wrong case](#pdb_id_case) | warning 5 | 5 |
| [Obsolete or superseded PDB IDs](#pdb_obsolete) | warning 34 | 34 |
| [DOI matches the PDB but the PubMed ID does not](#pubmed_id_conflict) | warning 1 | 1 |
| [PubMed ID points to a different paper than the citation title](#pubmed_title_mismatch) | warning 4 | 4 |
| [PDB IDs listed more than once in the same subgroup](#repeated_in_subgroup) | warning 5 | 5 |
| [Resolution differs from the PDB entry](#resolution_differs_from_pdb) | warning 65 | 65 |
| [Resolution text that is neither a number nor a method](#resolution_unparseable) | warning 2 | 2 |
| [Resolution conflicts with the PDB's experimental method](#resolution_vs_method) | warning 3 | 3 |
| [Entries listing themselves as related](#self_related) | warning 10 | 10 |
| [Unclosed <a> links](#unclosed_link) | warning 1 | 1 |
| [Taxonomic domain outside Eukaryota/Bacteria/Archaea/Viruses](#unknown_domain) | warning 6 | 6 |
| [Misspelled HTML entities](#unknown_entity) | warning 4 | 4 |
| [HTML entities missing ';'](#unterminated_entity) | warning 5 | 5 |
| [Unescaped '&' (makes the file invalid XML)](#bare_ampersand) | info 81 | 81 |
| [PDB IDs listed in more than one subgroup (usually intentional)](#cross_listing) | info 21 | 21 |
| [Subgroup names used more than once](#duplicate_subgroup_name) | info 2 | 2 |
| [Subgroups with no entries](#empty_subgroup) | info 1 | 1 |
| [Empty formatting tags](#empty_tag) | info 1 | 1 |
| [Zero-width characters](#invisible_characters) | info 2 | 2 |
| [PDB IDs not yet released](#pdb_unreleased) | info 6 | 6 |
| [Related entries not cited in the description](#related_not_described) | info 25 | 25 |
| [Resolution with several values](#resolution_multiple) | info 1 | 1 |

<a id="pdb_id_points_elsewhere"></a>
## PDB ID points to an unrelated structure

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| error | `G02-S007-P0016-D002` | 5LGO | PDB 5LGO is 'Trypsin inhibitors for the treatment of pancreatitis - cpd 15' (X-RAY DIFFRACTION), which does not match 'BamABCDE complex o… |
| error | `G03-S025-P0002-D001` | 8AT1 | PDB 8AT1 is 'CRYSTAL STRUCTURES OF ASPARTATE CARBAMOYLTRANSFERASE LIGATED WITH PHOSPHONOACETAMIDE, MALONATE, AND CTP OR ATP AT 2.8-ANGSTR… |
| error | `G03-S051-P0168` | 8HCD | PDB 8HCD is 'Crystal structure of mTREX1-DNA product complex (dNMP)' (X-RAY DIFFRACTION), which does not match 'endothelin ETB receptor -… |
| error | `G03-S094-P0001-D020` | 5OO6 | PDB 5OO6 is 'Complex of human nuclear cap-binding complex with ARS2 C-terminal peptide' (X-RAY DIFFRACTION), which does not match 'AcrABZ… |
| error | `G03-S098-P0007` | 6IJ6 | PDB 6IJ6 is 'Crystal structure of PETase S121E, D186H, R280A mutant from Ideonella sakaiensis' (X-RAY DIFFRACTION), which does not match … |
| error | `G03-S114-P0037` | 6SP7 | PDB 6SP7 is 'Crystal Structure of the VIM-2 Acquired Metallo-beta-Lactamase in Complex with Taniborbactam (VNRX-5133)' (X-RAY DIFFRACTION… |
| error | `G03-S146-P0011` | 7NK2 | PDB 7NK2 is '1918 H1N1 Viral influenza polymerase heterotrimer with Nb8202 core' (ELECTRON MICROSCOPY), which does not match 'Cytochrome … |

<a id="pdb_unknown"></a>
## PDB IDs the PDB does not recognise

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| error | `G03-S039-P0022-D004-R023` | 72LR | Not a current, unreleased or removed PDB entry |
| error | `G03-S132-P0005-D008-R001` | 8L39 | Not a current, unreleased or removed PDB entry |
| error | `G03-S132-P0005-D008-R002` | 8L29 | Not a current, unreleased or removed PDB entry |

<a id="bad_doi"></a>
## Malformed DOIs

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S029-P0006-D001` | 6XYN | primary DOI '103390/cryst10030149' |
| warning | `G03-S085-P0011` | 4AW6 | primary DOI 'doi:10.1126/science.1231513' |

<a id="broken_entity"></a>
## Malformed entity spellings

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S033-P0001-D001` | 6QM9 | description: malformed '& Aring;' (displayed as 'Å') |
| warning | `G03-S037-P0005` | 6RV2 | bibliography[0].authors: malformed '&uml;' (left as '¨'; probably ü) |
| warning | `G03-S094-P0001-D011` | 3AOB | description: malformed 'Aring;' (displayed as 'Å') |

<a id="citation_differs_from_pdb"></a>
## Primary citation differs from the PDB entry's citation

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G01-S001-P0001` | 1PRH | Primary citation 'Crystal structure of an integral membrane light-harvesting complex from photosynthetic bacteria.' differs from the PDB'… |
| warning | `G01-S001-P0001-D002` | 1CQE | Primary citation 'The 3.1 A X-ray crystal structure of the integral membrane enzyme prostaglandin H2 synthase-1.' differs from the PDB's … |
| warning | `G02-S005-P0074` | 8YPU | Primary citation 'Gastric mucous neck cell and intestinal goblet cell phenotypes in gastric adenocarcinoma.' differs from the PDB's 'A hi… |
| warning | `G02-S007-P0001` | 4QKY | Primary citation 'Structure of the membrane protein FhaC: a member of the Omp85-TpsB transporter superfamily.' differs from the PDB's 'Co… |
| warning | `G03-S027-P0017` | 6CSM | Primary citation 'Crystal structure of the natural anion-conducting channelrhodopsin GtACR1.' differs from the PDB's 'Structural mechanis… |
| warning | `G03-S033-P0002-D006` | 8XLR | Primary citation 'Structural mechanism of voltage-gated sodium channel slow inactivation.' differs from the PDB's 'Tamsulosin ameliorates… |
| warning | `G03-S036-P0014` | 8GRO | Primary citation 'Mechanism underlying delayed rectifying in human voltage-mediated activation Eag2 channel.' differs from the PDB's 'A m… |
| warning | `G03-S039-P0012-D004` | 9Y2A | Primary citation 'Shape-shifting conotoxins reveal divergent pore-targeting mechanisms in nicotinic receptors.' differs from the PDB's 'N… |
| warning | `G03-S039-P0024-D001` | 9P6B | Primary citation 'Structural basis for excitatory neuropeptide signaling.' differs from the PDB's 'Structural basis of the inhibition of … |
| warning | `G03-S039-P0034-D002` | 6DO7 | Primary citation 'Structural bases of TRP channel TRPV6 allosteric modulation by 2-APB.' differs from the PDB's 'A ligand-induced structu… |
| warning | `G03-S045-P0007` | 9J7D | Primary citation 'Molecular basis of pyruvate transport and inhibition of the human mitochondrial pyruvate carrier.' differs from the PDB… |
| warning | `G03-S048-P0013` | 2MAW | Primary citation 'NMR structures of the human α7 nAChR transmembrane domain and associated anesthetic binding sites.' differs from the PD… |
| warning | `G03-S048-P0018` | 8EIS | Primary citation 'Sensory specializations drive octopus and squid behaviour.' differs from the PDB's 'Structural basis of sensory recepto… |
| warning | `G03-S048-P0022-D008` | 43SN | Primary citation 'Cryo-EM structure of the bicarbonate receptor GPR30.' differs from the PDB's 'Screening, synthesis and structural chara… |
| warning | `G03-S051-P0024-D001` | 22GK | Primary citation 'Molecular mechanisms of naturally encoded signaling bias at the complement anaphylatoxin receptors.' differs from the P… |
| warning | `G03-S051-P0025-D004` | 22GL | Primary citation 'Molecular mechanisms of naturally encoded signaling bias at the complement anaphylatoxin receptors.' differs from the P… |
| warning | `G03-S085-P0001-D010` | 2LEP | Primary citation 'Insights into the effect of detergents on the full-length rhomboid protease from Pseudomonas aeruginosa and its cytosol… |
| warning | `G03-S094-P0019` | 4MT1 | Primary citation 'Crystal structure of the Neisseria gonorrhoeae MtrD inner membrane multidrug efflux pump.' differs from the PDB's 'Crys… |
| warning | `G03-S124-P0001` | 5U73 | Primary citation 'Structure of human Niemann-Pick C1 protein.' differs from the PDB's '3.3 angstrom structure of Niemann-Pick C1 protein … |
| warning | `G03-S149-P0001` | 2PPS | Primary citation 'Photosystem I of Synechococcus elongatus at 4 Å resolution: Comprehensive structure analysis.' differs from the PDB's '… |
| warning | `G03-S153-P0001` | 1PRC | Primary citation 'Structure of the protein subunits in the photosynthetic reaction centre of Rhodospeudomonas viridis at 3 Å resolution.'… |
| warning | `G03-S153-P0001-D003` | 5O4C | Primary citation 'Ultrafast structural changes within a photosynthetic reaction centre.' differs from the PDB's 'From Macrocrystals to Mi… |
| info | `G03-S146-P0002` | 1AR1 | PDB primary citation is listed here as a secondary reference |
| info | `G03-S153-P0002` | 4RCR | PDB primary citation is listed here as a secondary reference |

<a id="citation_ids_wrong"></a>
## Citation title matches the PDB's, but its PubMed ID/DOI do not

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S027-P0001-D001` | 1AT9 | Same title as the PDB's citation, but the PubMed ID/DOI differ (mpstruc 8676377 / —, PDB 9296502 / 10.1038/38323) |
| warning | `G03-S133-P0001-D007` | 2BY4 | Same title as the PDB's citation, but the PubMed ID/DOI differ (mpstruc 16150713 / —, PDB 16412648 / 10.1016/j.bmc.2005.12.001) |

<a id="citation_preprint"></a>
## mpstruc cites a preprint that has since been published

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S051-P0098-D001` | 8UGY | mpstruc cites the preprint; the published version is Sci Adv 2024 (DOI 10.1126/sciadv.adk4855, PMID 38630816) |
| warning | `G03-S054-P0007-D005` | 8T7H | mpstruc cites the preprint; the published version is Nature 2024 (DOI 10.1038/s41586-024-07327-x, PMID 38632403) |

<a id="described_not_related"></a>
## Description cites a structure missing from the related list

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G02-S003-P0026` | 3VZT | Description cites 3ZVU but its related list does not include it |
| warning | `G03-S037-P0112` | 8ON8 | Description cites 8ON7 but its related list does not include it (listed elsewhere: 8ON7 at G03-S039-P0024-D001-R001) |
| warning | `G03-S041-P0001-D020` | 9B5Z | Description cites 8QOJ, 8R7Q, 8R7R but its related list does not include them (listed elsewhere: 8QOJ at G03-S040-P0002-D002-R001; 8R7Q a… |
| warning | `G03-S041-P0047-D002` | 9EGT | Description cites 9EGS but its related list does not include it |
| warning | `G03-S051-P0023` | 8I97 | Description cites 8J6D but its related list does not include it |
| warning | `G03-S051-P0027` | 9WDI | Description cites 25IF but its related list does not include it |
| warning | `G03-S051-P0094-D002` | 7W0L | Description cites 7W0N but its related list does not include it |
| warning | `G03-S051-P0116` | 6RZ6 | Description cites 6R7Z but its related list does not include it (listed elsewhere: 6R7Z at G03-S033-P0005-R003) |
| warning | `G03-S051-P0146-D001` | 7F8U | Description cites 7F8W but its related list does not include it (listed elsewhere: 7F8W at G03-S051-P0147-R001) |
| warning | `G03-S061-P0001-D002` | 9GS2 | Description cites 9GSN but its related list does not include it (listed elsewhere: 9GSN at G03-S061-P0001-D001-R001) |
| warning | `G03-S102-P0042-D001` | 6LZ0 | Description cites 7CKR, 7CKO, 7DA5 but its related list does not include them (listed elsewhere: 7CKR at G03-S102-P0042-R001; 7CKO at G03… |
| warning | `G03-S102-P0077` | 8J74 | Description cites 8J77 but its related list does not include it |
| warning | `G03-S131-P0008` | 9SIE | Description cites 9XB9 but its related list does not include it (listed elsewhere: 9XB9 at G03-S037-P0098-D003-R003) |
| warning | `G03-S139-P0001` | 2VPZ | Description cites 2PVY, 2PVW, 2PVX but its related list does not include them |

<a id="doi_conflict"></a>
## PubMed ID matches the PDB but the DOI does not

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G01-S002-P0003` | 1LNH | PubMed ID matches the PDB, but DOI '10.1002/(sici)1097-0134(199709)29:1' differs from the PDB's '10.1002/(sici)1097-0134(199709)29:1<15::… |
| warning | `G02-S003-P0006-D003` | 1GFM | PubMed ID matches the PDB, but DOI '10.1074/jbc.271.34.20669' differs from the PDB's '10.1074/jbc.271.34.20676' |
| warning | `G03-S037-P0003` | 3UKM | PubMed ID matches the PDB, but DOI '10.1126/science.121327' differs from the PDB's '10.1126/science.1213274' (truncated?); PubMed has '10… |
| info | `G03-S051-P0138` | 8IOC | The PDB's DOI '10.1107/s0907444909042073' looks wrong; mpstruc's '10.1038/s41421-023-00586-4' agrees with PubMed |
| info | `G03-S051-P0140` | 8INR | The PDB's DOI '10.1107/s0907444909042073' looks wrong; mpstruc's '10.1038/s41421-023-00586-4' agrees with PubMed |
| info | `G03-S069-P0002` | 8RF0 | The PDB's DOI '10.1038/s41467-020-17921-y' looks wrong; mpstruc's '10.1038/s41467-024-45415-8' agrees with PubMed |
| info | `G03-S094-P0001-D030` | 9BFH | The PDB's DOI '10.1073/pnas.1901346116' looks wrong; mpstruc's '10.1128/mbio.02492-23' agrees with PubMed |
| info | `G03-S149-P0014-D011` | 11YO | The PDB's DOI '10.1073/pnas.202398211' looks wrong; mpstruc's '10.1073/pnas.2023982118' agrees with PubMed |

<a id="malformed_markup"></a>
## Unexpected or malformed tags

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S038-P0017-D003` | 8TK8 | description: unexpected tag(s) <ahref> |
| warning | `G03-S038-P0017-D004` | 9YKK | description: unexpected tag(s) <super> |

<a id="missing_field"></a>
## Empty name, species or taxonomic domain

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G02-S003-P0032` | 5LDV | taxonomicDomain is empty |
| warning | `G02-S005-P0042` | 6BPN | species is empty |
| warning | `G03-S037-P0001-D017` | 6W0A | taxonomicDomain is empty |
| warning | `G03-S045-P0004-D002` | 9W4L | species is empty |
| warning | `G03-S051-P0062-D001` | 6IQL | species is empty |
| warning | `G03-S051-P0095-D004` | 9PJK | taxonomicDomain is empty |
| warning | `G03-S051-P0183` | 8XOG | taxonomicDomain is empty |
| warning | `G03-S051-P0209` | 8W1Z | species is empty |
| warning | `G03-S102-P0077` | 8J74 | taxonomicDomain is empty |
| warning | `G03-S106-P0001-D001` | 5XLS | taxonomicDomain is empty |
| warning | `G03-S112-P0003-D004` | 9QZQ | taxonomicDomain is empty |
| warning | `G03-S115-P0003` | 9PYM | species is empty |
| warning | `G03-S115-P0003` | 9PYM | taxonomicDomain is empty |
| warning | `G03-S130-P0007` | 5I1M | taxonomicDomain is empty |

<a id="pdb_id_case"></a>
## PDB IDs in the wrong case

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G02-S007-P0016-D010-R001` | 8ADi | '8ADi' should be '8ADI' (PDB IDs are upper case) |
| warning | `G03-S038-P0015-D001-R001` | 7vMS | '7vMS' should be '7VMS' (PDB IDs are upper case) |
| warning | `G03-S039-P0014-D001-R001` | 8Si4 | '8Si4' should be '8SI4' (PDB IDs are upper case) |
| warning | `G03-S039-P0014-D001-R004` | 8Si3 | '8Si3' should be '8SI3' (PDB IDs are upper case) |
| warning | `G03-S094-P0030-R002` | 7RPj | '7RPj' should be '7RPJ' (PDB IDs are upper case) |

<a id="pdb_obsolete"></a>
## Obsolete or superseded PDB IDs

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G01-S004-P0001-D001-R003` | 1OJB | superseded: 1OJB → 2XFU |
| warning | `G01-S005-P0003-R001` | 4GGI | superseded: 4GGI → 4J6E |
| warning | `G01-S015-P0001` | 5W7L | superseded: 5W7L → 8G1N |
| warning | `G03-S003-P0001` | 5G1J | obsolete: 5G1J (no replacement) |
| warning | `G03-S023-P0001` | 3WXV | superseded: 3WXV → 6KS0 |
| warning | `G03-S023-P0001-R001` | 3WXW | superseded: 3WXW → 6KS1 |
| warning | `G03-S027-P0024-D001-R009` | 6REW | superseded: 6REW → 6YC3 |
| warning | `G03-S031-P0001-R001` | 3HD9 | superseded: 3HD9 → 3IPD |
| warning | `G03-S038-P0013-D001` | 3J8E | superseded: 3J8E → 5TB0 |
| warning | `G03-S041-P0025-D001` | 3HGC | superseded: 3HGC → 4NYK |
| warning | `G03-S042-P0001-R002` | 5A41 | superseded: 5A41 → 5FXB |
| warning | `G03-S051-P0147-D001-R001` | 7XOX | superseded: 7XOX → 8IA7 |
| warning | `G03-S060-P0002-D002-R001` | 4FTP | superseded: 4FTP → 6LSC |
| warning | `G03-S085-P0013-D002` | 4UPC | superseded: 4UPC → 5A63 |
| warning | `G03-S094-P0010-R003` | 3VVQ | superseded: 3VVQ → 3WBN |
| warning | `G03-S101-P0003-D003-R003` | 4OYG | superseded: 4OYG → 5CFY |
| warning | `G03-S105-P0001-R001` | 2JLO | superseded: 2JLO → 4D1B |
| warning | `G03-S113-P0008` | 6FFV | superseded: 6FFV → 8C7P |
| warning | `G03-S114-P0036` | 6AN7 | superseded: 6AN7 → 6OIH |
| warning | `G03-S114-P0043` | 5TSI | superseded: 5TSI → 5UAR |
| warning | `G03-S116-P0015` | 4J05 | superseded: 4J05 → 7SP5 |
| warning | `G03-S133-P0001-D003-R001` | 2ZDB | superseded: 2ZDB → 3B02 |
| warning | `G03-S133-P0001-D016-R001` | 2EAS | superseded: 2EAS → 4YCL |
| warning | `G03-S133-P0006-D006-R001` | 7DDG | superseded: 7DDG → 7EVX |
| warning | `G03-S133-P0006-D006-R005` | 7DDJ | superseded: 7DDJ → 7WYT |
| warning | `G03-S133-P0013` | 3B8C | superseded: 3B8C → 5KSD |
| warning | `G03-S140-P0015` | 4P6V | superseded: 4P6V → 8ACY |
| warning | `G03-S143-P0001` | 1FUM | superseded: 1FUM → 1L0V |
| warning | `G03-S145-P0001-R001` | 1UM3 | superseded: 1UM3 → 1VF5 |
| warning | `G03-S149-P0014-D003` | 3BZ1 | superseded: 3BZ1 → 4V62 |
| warning | `G03-S149-P0014-D003-R001` | 3BZ2 | superseded: 3BZ2 → 4V62 |
| warning | `G03-S149-P0014-D004-R001` | 4Q54 | superseded: 4Q54 → 4RVY |
| warning | `G03-S149-P0015-D002` | 3ARC | superseded: 3ARC → 3WU2 |
| warning | `G03-S151-P0010` | 3CJU | superseded: 3CJU → 3EGV |

<a id="pubmed_id_conflict"></a>
## DOI matches the PDB but the PubMed ID does not

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S142-P0001-D001` | 3RKO | DOI matches the PDB, but PubMed ID 21831629 differs from the PDB's 21822288 (Respiratory complex I: 'steam engine' of the cell?) |

<a id="pubmed_title_mismatch"></a>
## PubMed ID points to a different paper than the citation title

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S027-P0001-D001` | 1AT9 | Primary PubMed ID 8676377 is 'Electron-crystallographic refinement of the structure of bacteriorhodopsin.' (J Mol Biol 1996), but the cit… |
| warning | `G03-S051-P0035-D005` | 3P0G | Primary PubMed ID 21228869 is 'Structure of a nanobody-stabilized active state of the β(2) adrenoceptor.' (Nature 2011), but the citation… |
| warning | `G03-S133-P0001-D007` | 2BY4 | Primary PubMed ID 16150713 is 'Structural role of countertransport revealed in Ca(2+) pump crystal structure in the absence of Ca(2+).' (… |
| warning | `G03-S142-P0001-D001` | 3RKO | Primary PubMed ID 21831629 is 'Respiratory complex I: 'steam engine' of the cell?' (Curr Opin Struct Biol 2011), but the citation title i… |

<a id="repeated_in_subgroup"></a>
## PDB IDs listed more than once in the same subgroup

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S039-P0034-D003` | 7D2K | Listed 2 times within subgroup G03-S039: G03-S039-P0034-D003, G03-S039-P0035-D003-R003 |
| warning | `G03-S051-P0024-D001-R001` | 22GL | Listed 2 times within subgroup G03-S051: G03-S051-P0024-D001-R001, G03-S051-P0025-D004 |
| warning | `G03-S051-P0146-D001-R002` | 7F8V | Listed 2 times within subgroup G03-S051: G03-S051-P0146-D001-R002, G03-S051-P0147 |
| warning | `G03-S132-P0001` | 1BMF | Listed 2 times within subgroup G03-S132: G03-S132-P0001, G03-S132-P0005-R002 |
| warning | `G03-S132-P0003` | 1A91 | Listed 2 times within subgroup G03-S132: G03-S132-P0003, G03-S132-P0005-R003 |

<a id="resolution_differs_from_pdb"></a>
## Resolution differs from the PDB entry

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G01-S002-P0004-D001` | 3FG4 | mpstruc resolution 1.85 Å, PDB reports 2.31 Å |
| warning | `G02-S001-P0015` | 5FMW | mpstruc resolution 8.0 Å, PDB reports 6.7 Å |
| warning | `G02-S003-P0006-D013` | 6WTZ | mpstruc resolution 2.54 Å, PDB reports 3.15 Å |
| warning | `G02-S003-P0028` | 4FQE | mpstruc resolution 2.1 Å, PDB reports 1.93 Å |
| warning | `G02-S004-P0007` | 3T0S | mpstruc resolution 2.1 Å, PDB reports 2.2 Å |
| warning | `G02-S005-P0006` | 5AZO | mpstruc resolution 1.69 Å, PDB reports 2.7 Å |
| warning | `G02-S005-P0063` | 6GIE | mpstruc resolution 2.2 Å, PDB reports 2.1 Å |
| warning | `G02-S007-P0001` | 4QKY | mpstruc resolution 3.15 Å, PDB reports 2.9 Å |
| warning | `G02-S007-P0006` | 2QCZ | mpstruc resolution 2.2 Å, PDB reports 2.7 Å |
| warning | `G02-S007-P0006-D006` | 6FSU | mpstruc resolution 2.2 Å, PDB reports 2.6 Å |
| warning | `G02-S007-P0016-D002` | 5LGO | mpstruc resolution 4.9 Å, PDB reports 1.12 Å |
| warning | `G02-S008-P0001-D001` | 2JK4 | mpstruc resolution 4.0 Å, PDB reports 4.1 Å |
| warning | `G03-S027-P0001-D001` | 1AT9 | mpstruc resolution 3.0 Å, PDB reports 2.8, 2.8 Å |
| warning | `G03-S030-P0001-D004` | 9CQ4 | mpstruc resolution 3.7 Å, PDB reports 3.27 Å |
| warning | `G03-S032-P0003` | 4P79 | mpstruc resolution 2.61 Å, PDB reports 2.4 Å |
| warning | `G03-S035-P0002` | 3MP7 | mpstruc resolution 3.1 Å, PDB reports 2.9 Å |
| warning | `G03-S037-P0007` | 3UM7 | mpstruc resolution 3.8 Å, PDB reports 3.31 Å |
| warning | `G03-S037-P0016-D001` | 8TEO | mpstruc resolution 3.0 Å, PDB reports 2.39 Å |
| warning | `G03-S037-P0044` | 3SYO | mpstruc resolution 3.6 Å, PDB reports 3.54 Å |
| warning | `G03-S037-P0055-D003` | 6C1E | mpstruc resolution 2.52 Å, PDB reports 2.86 Å |
| warning | `G03-S037-P0094` | 8INZ | mpstruc resolution 3.19 Å, PDB reports 2.72 Å |
| warning | `G03-S038-P0014-D002` | 6JG3 | mpstruc resolution 6.3 Å, PDB reports 6.1 Å |
| warning | `G03-S038-P0017-D003` | 8TK8 | mpstruc resolution 3.49 Å, PDB reports 2.7 Å |
| warning | `G03-S039-P0025` | 5AN8 | mpstruc resolution 4.0 Å, PDB reports 3.8 Å |
| warning | `G03-S040-P0002` | 7XKK | mpstruc resolution 2.2 Å, PDB reports 3.2 Å |
| warning | `G03-S041-P0001-D020` | 9B5Z | mpstruc resolution 2.53 Å, PDB reports 2.71 Å |
| warning | `G03-S045-P0006` | 8XDG | mpstruc resolution 2.4 Å, PDB reports 3.1 Å |
| warning | `G03-S048-P0002-D005` | 8QQM | mpstruc resolution 2.9 Å, PDB reports 4.7 Å |
| warning | `G03-S048-P0022-D008` | 43SN | mpstruc resolution 3.15 Å, PDB reports 2.71 Å |
| warning | `G03-S051-P0038-D030` | 9T9P | mpstruc resolution 2.6 Å, PDB reports 3.0 Å |
| warning | `G03-S051-P0058-D003` | 7T6V | mpstruc resolution 2.9 Å, PDB reports 3.1 Å |
| warning | `G03-S051-P0058-D004` | 8ZBW | mpstruc resolution 3.17 Å, PDB reports 2.58 Å |
| warning | `G03-S051-P0121-D002` | 9XFI | mpstruc resolution 2.64 Å, PDB reports 2.99 Å |
| warning | `G03-S051-P0168` | 8HCD | mpstruc resolution 2.99 Å, PDB reports 2.0 Å |
| warning | `G03-S051-P0191-D001` | 8ZX4 | mpstruc resolution 3.03 Å, PDB reports 2.85 Å |
| warning | `G03-S051-P0192-D003` | 9M88 | mpstruc resolution 3.83 Å, PDB reports 3.1 Å |
| warning | `G03-S053-P0005-D002` | 8JMT | mpstruc resolution 3.3 Å, PDB reports 3.36 Å |
| warning | `G03-S055-P0001` | 7AD3 | mpstruc resolution 3.5 Å, PDB reports 3.3 Å |
| warning | `G03-S060-P0002-D002` | 4ENE | mpstruc resolution 2.5 Å, PDB reports 2.4 Å |
| warning | `G03-S062-P0002` | 9WE7 | mpstruc resolution 3.27 Å, PDB reports 3.05 Å |
| warning | `G03-S067-P0002` | 9I5L | mpstruc resolution 3.2 Å, PDB reports 3.4 Å |
| warning | `G03-S085-P0013` | 4UIS | mpstruc resolution 4.32 Å, PDB reports 4.4 Å |
| warning | `G03-S094-P0001-D020` | 5OO6 | mpstruc resolution 5.9 Å, PDB reports 2.8 Å |
| warning | `G03-S097-P0004` | 2Q7M | mpstruc resolution 4.0 Å, PDB reports 4.25 Å |
| warning | `G03-S098-P0007` | 6IJ6 | mpstruc resolution 2.23 Å, PDB reports 1.95 Å |
| warning | `G03-S101-P0004-D003` | 8AFA | mpstruc resolution 3.16 Å, PDB reports 3.25 Å |
| warning | `G03-S102-P0060` | 6XPD | mpstruc resolution 4.1 Å, PDB reports 3.8 Å |
| warning | `G03-S102-P0077-D001` | 8ZQR | mpstruc resolution 2.6 Å, PDB reports 3.6 Å |
| warning | `G03-S102-P0091-D001` | 9IJZ | mpstruc resolution 3.79 Å, PDB reports 2.91 Å |
| warning | `G03-S113-P0005-D001` | 7NNU | mpstruc resolution 2.2 Å, PDB reports 2.7 Å |
| warning | `G03-S114-P0012-D001` | 6QEX | mpstruc resolution 3.5 Å, PDB reports 3.6 Å |
| warning | `G03-S114-P0024-D001` | 7Y48 | mpstruc resolution 2.58 Å, PDB reports 2.85 Å |
| warning | `G03-S114-P0037` | 6SP7 | mpstruc resolution 3.2 Å, PDB reports 1.8 Å |
| warning | `G03-S119-P0003` | 5D91 | mpstruc resolution 3.6 Å, PDB reports 2.501 Å |
| warning | `G03-S132-P0012-D001` | 8AP8 | mpstruc resolution 3.2 Å, PDB reports 3.7 Å |
| warning | `G03-S133-P0001-D006` | 2AGV | mpstruc resolution 3.0 Å, PDB reports 2.4 Å |
| warning | `G03-S133-P0004-D002` | 7W7T | mpstruc resolution 2.9 Å, PDB reports 3.4 Å |
| warning | `G03-S140-P0006-D001` | 2ZUQ | mpstruc resolution 3.4 Å, PDB reports 3.3 Å |
| warning | `G03-S140-P0017-D001` | 8UCD | mpstruc resolution 2.0 Å, PDB reports 3.0 Å |
| warning | `G03-S142-P0007-D004` | 7PSA | mpstruc resolution 3.04 Å, PDB reports 3.4 Å |
| warning | `G03-S142-P0009` | 4WZ7 | mpstruc resolution 3.8 Å, PDB reports 3.6 Å |
| warning | `G03-S146-P0001-D007` | 6NMP | mpstruc resolution 1.75 Å, PDB reports 2.9 Å |
| warning | `G03-S146-P0005-D001` | 9VDJ | mpstruc resolution 2.36 Å, PDB reports 2.3 Å |
| warning | `G03-S146-P0007-D001` | 6WTI | mpstruc resolution 2.2 Å, PDB reports 2.38 Å |
| warning | `G03-S146-P0011` | 7NK2 | mpstruc resolution 2.5 Å, PDB reports 4.84 Å |

<a id="resolution_unparseable"></a>
## Resolution text that is neither a number nor a method

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S030-P0011` | 8J3V | resolution '-' is neither a number nor a method |
| warning | `G03-S037-P0018-D004` | 8SD3 | resolution 'HEK293 cells' is neither a number nor a method |

<a id="resolution_vs_method"></a>
## Resolution conflicts with the PDB's experimental method

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S025-P0002-D001` | 8AT1 | mpstruc says 'NMR structure', but the PDB entry is X-RAY DIFFRACTION at 2.8 Å |
| warning | `G03-S039-P0034-D002` | 6DO7 | mpstruc gives 3.45 Å, but the PDB entry is SOLUTION NMR (no resolution) |
| warning | `G03-S140-P0006-D004` | 2LTQ | mpstruc gives 1.35 Å, but the PDB entry is SOLID-STATE NMR (no resolution) |

<a id="self_related"></a>
## Entries listing themselves as related

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G02-S008-P0001-D002` | 6TIQ | Entry lists itself as related |
| warning | `G03-S027-P0024-D003` | 4XTL | Entry lists itself as related |
| warning | `G03-S051-P0196` | 8XFT | Entry lists itself as related |
| warning | `G03-S077-P0001` | 6BML | Entry lists itself as related |
| warning | `G03-S082-P0003` | 9L0I | Entry lists itself as related |
| warning | `G03-S102-P0060` | 6XPD | Entry lists itself as related |
| warning | `G03-S106-P0001-D002` | 8OMZ | Entry lists itself as related |
| warning | `G03-S109-P0016-D001` | 9GIV | Entry lists itself as related |
| warning | `G03-S132-P0005-D004` | 6B2Z | Entry lists itself as related |
| warning | `G03-S154-P0002` | 7UEB | Entry lists itself as related |

<a id="unclosed_link"></a>
## Unclosed <a> links

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G03-S039-P0004-D001` | 7B0S | description: 1 <a> vs 0 </a> |

<a id="unknown_domain"></a>
## Taxonomic domain outside Eukaryota/Bacteria/Archaea/Viruses

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G02-S002-P0001` | 6X1K | taxonomicDomain 'Unclassified' |
| warning | `G03-S002-P0001` | 4P6J | taxonomicDomain 'Unclassified' |
| warning | `G03-S002-P0002` | 6B87 | taxonomicDomain 'Unclassified' |
| warning | `G03-S002-P0003` | 6MQU | taxonomicDomain 'Unclassified' |
| warning | `G03-S002-P0004` | 6TJ1 | taxonomicDomain 'Unclassified' |
| warning | `G03-S027-P0015` | 6CSN | taxonomicDomain 'Unclassified' |

<a id="unknown_entity"></a>
## Misspelled HTML entities

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G02-S003-P0029` | 4V3G | description: &alfa; is not an HTML entity (displayed as 'α') |
| warning | `G03-S038-P0001-D001` | 5KLB | description: &Ariing; is not an HTML entity (displayed as 'Å') |
| warning | `G03-S133-P0021-D002` | 7KY9 | description: &Ariing; is not an HTML entity (displayed as 'Å') |
| warning | `G03-S140-P0019` | 6WXR | description: &Ariing; is not an HTML entity (displayed as 'Å') |

<a id="unterminated_entity"></a>
## HTML entities missing ';'

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| warning | `G02-S001-P0001-D003` | 6U49 | description: '&Aring' is missing its ';' |
| warning | `G03-S027-P0001-D018` | 6G7H | description: '&Aring' is missing its ';' |
| warning | `G03-S054-P0007-D003` | 6FFI | bibliography[0].authors: '&eacute' is missing its ';' |
| warning | `G03-S142-P0017` | 7ARB | bibliography[0].authors: '&Ouml' is missing its ';' |
| warning | `G03-S142-P0018` | 7ARD | bibliography[0].authors: '&Ouml' is missing its ';' |

<a id="bare_ampersand"></a>
## Unescaped '&' (makes the file invalid XML)

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| info | `G01-S005-P0002` | 5B49 | name: unescaped '&' (invalid XML) |
| info | `G01-S013-P0006` | 6EJI | description: unescaped '&' (invalid XML) |
| info | `G01-S030-P0003` | 6Z7A | name: unescaped '&' (invalid XML) |
| info | `G02-S007-P0005` | 4K3C | name: unescaped '&' (invalid XML) |
| info | `G02-S008-P0007` | 6WUH | name: unescaped '&' (invalid XML) |
| info | `G02-S008-P0007-D001` | 9NK6 | name: unescaped '&' (invalid XML) |
| info | `G02-S011-P0003` | 6HCG | description: unescaped '&' (invalid XML) |
| info | `G03-S006-P0001-D001` | 9NM2 | name: unescaped '&' (invalid XML) |
| info | `G03-S027-P0014-D002` | 7C86 | name: unescaped '&' (invalid XML) |
| info | `G03-S033-P0001-D001` | 6QM9 | description: unescaped '&' (invalid XML) |
| info | `G03-S033-P0002-D005` | 8QZC | name: unescaped '&' (invalid XML) |
| info | `G03-S033-P0004` | 6QPC | name: unescaped '&' (invalid XML) |
| info | `G03-S033-P0004-D001` | 6P46 | name: unescaped '&' (invalid XML) |
| info | `G03-S033-P0004-D002` | 8B8G | name: unescaped '&' (invalid XML) |
| info | `G03-S033-P0004-D003` | 8TAG | name: unescaped '&' (invalid XML) |
| info | `G03-S033-P0005` | 6R7X | name: unescaped '&' (invalid XML) |
| info | `G03-S037-P0001-D003` | 2BOB | name: unescaped '&' (invalid XML) |
| info | `G03-S037-P0001-D013` | 2JK5 | name: unescaped '&' (invalid XML) |
| info | `G03-S037-P0001-D017` | 6W0A | description: unescaped '&' (invalid XML) |
| info | `G03-S037-P0009` | 6CQ6 | description: unescaped '&' (invalid XML) |
| info | `G03-S037-P0011-D003` | 6E1M | description: unescaped '&' (invalid XML) |
| info | `G03-S037-P0022` | 7W3Y | description: unescaped '&' (invalid XML) |
| info | `G03-S037-P0107` | 8JET | description: unescaped '&' (invalid XML) |
| info | `G03-S037-P0107-D001` | 8WTZ | description: unescaped '&' (invalid XML) |
| info | `G03-S038-P0013-D003` | 6FOO | name: unescaped '&' (invalid XML) |
| info | `G03-S038-P0015-D001` | 7VML | name: unescaped '&' (invalid XML) |
| info | `G03-S038-P0016-D001` | 7LHE | description: unescaped '&' (invalid XML) |
| info | `G03-S039-P0016` | 6BPQ | name: unescaped '&' (invalid XML) |
| info | `G03-S039-P0016-D001` | 6NR3 | name: unescaped '&' (invalid XML) |
| info | `G03-S039-P0016-D002` | 6O6A | name: unescaped '&' (invalid XML) |
| info | `G03-S039-P0016-D003` | 8E4Q | name: unescaped '&' (invalid XML) |
| info | `G03-S039-P0017` | 8E4P | name: unescaped '&' (invalid XML) |
| info | `G03-S039-P0017-D001` | 9B6D | name: unescaped '&' (invalid XML) |
| info | `G03-S039-P0018` | 8BDC | name: unescaped '&' (invalid XML) |
| info | `G03-S039-P0039` | 5TJA | bibliography[0].authors: unescaped '&' (invalid XML) |
| info | `G03-S039-P0041` | 9B6I | name: unescaped '&' (invalid XML) |
| info | `G03-S041-P0001-D006` | 5KBS | species: unescaped '&' (invalid XML) |
| info | `G03-S041-P0001-D017` | 7TNJ | name: unescaped '&' (invalid XML) |
| info | `G03-S041-P0006-D002` | 5UOW | name: unescaped '&' (invalid XML) |
| info | `G03-S041-P0071` | 7LZH | description: unescaped '&' (invalid XML) |
| info | `G03-S048-P0014-D005` | 6HUG | description: unescaped '&' (invalid XML) |
| info | `G03-S051-P0001-D022` | 7MT9 | description: unescaped '&' (invalid XML) |
| info | `G03-S051-P0033-D005` | 6H7J | name: unescaped '&' (invalid XML) |
| info | `G03-S051-P0033-D005` | 6H7J | description: unescaped '&' (invalid XML) |
| info | `G03-S051-P0035-D010` | 6E67 | description: unescaped '&' (invalid XML) |
| info | `G03-S051-P0037` | 5UEN | description: unescaped '&' (invalid XML) |
| info | `G03-S051-P0049` | 7XBW | description: unescaped '&' (invalid XML) |
| info | `G03-S051-P0050` | 4MBS | description: unescaped '&' (invalid XML) |
| info | `G03-S051-P0050-D002` | 7F1Q | name: unescaped '&' (invalid XML) |
| info | `G03-S051-P0072-D001` | 6OIJ | expressedInSpecies: unescaped '&' (invalid XML) |
| info | `G03-S051-P0073-D002` | 6OIK | expressedInSpecies: unescaped '&' (invalid XML) |
| info | `G03-S051-P0073-D003` | 6U1N | species: unescaped '&' (invalid XML) |
| info | `G03-S051-P0083-D005` | 5XRA | description: unescaped '&' (invalid XML) |
| info | `G03-S051-P0145` | 7JHJ | species: unescaped '&' (invalid XML) |
| info | `G03-S054-P0010-D005` | 7M3F | description: unescaped '&' (invalid XML) |
| info | `G03-S054-P0011` | 7DD5 | description: unescaped '&' (invalid XML) |
| info | `G03-S059-P0002-D001` | 6NT9 | species: unescaped '&' (invalid XML) |
| info | `G03-S062-P0001` | 7VHP | name: unescaped '&' (invalid XML) |
| info | `G03-S074-P0001-D001` | 5XHQ | bibliography[0].authors: unescaped '&' (invalid XML) |
| info | `G03-S088-P0003` | 5N77 | description: unescaped '&' (invalid XML) |
| info | `G03-S094-P0026-D001` | 6OEU | description: unescaped '&' (invalid XML) |
| info | `G03-S094-P0026-D003` | 6RVD | description: unescaped '&' (invalid XML) |
| info | `G03-S094-P0027-D002` | 7V6Y | description: unescaped '&' (invalid XML) |
| info | `G03-S101-P0001-D006` | 3QS4 | description: unescaped '&' (invalid XML) |
| info | `G03-S101-P0003-D008` | 6UWF | name: unescaped '&' (invalid XML) |
| info | `G03-S102-P0019` | 6YUP | name: unescaped '&' (invalid XML) |
| info | `G03-S102-P0066-D002` | 6VRH | name: unescaped '&' (invalid XML) |
| info | `G03-S111-P0010` | 7D3U | name: unescaped '&' (invalid XML) |
| info | `G03-S111-P0011` | 7QRU | name: unescaped '&' (invalid XML) |
| info | `G03-S114-P0007-D004` | 8DMO | name: unescaped '&' (invalid XML) |
| info | `G03-S114-P0007-D006` | 8TSO | name: unescaped '&' (invalid XML) |
| info | `G03-S114-P0007-D006` | 8TSO | description: unescaped '&' (invalid XML) |
| info | `G03-S114-P0011` | 6FN1 | species: unescaped '&' (invalid XML) |
| info | `G03-S114-P0044-D003` | 7SVD | name: unescaped '&' (invalid XML) |
| info | `G03-S114-P0044-D003` | 7SVD | description: unescaped '&' (invalid XML) |
| info | `G03-S114-P0044-D009` | 9DW5 | description: unescaped '&' (invalid XML) |
| info | `G03-S130-P0001-D008` | 6PE4 | species: unescaped '&' (invalid XML) |
| info | `G03-S130-P0002` | 6VQ6 | name: unescaped '&' (invalid XML) |
| info | `G03-S130-P0002` | 6VQ6 | species: unescaped '&' (invalid XML) |
| info | `G03-S141-P0003` | 7S9Y | name: unescaped '&' (invalid XML) |
| info | `G03-S151-P0002` | 5ZJI | name: unescaped '&' (invalid XML) |

<a id="cross_listing"></a>
## PDB IDs listed in more than one subgroup (usually intentional)

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| info | `G01-S010-P0001` | 3HYW | Listed in 2 subgroups (noted in description) |
| info | `G01-S010-P0001-R001` | 3HYX | Listed in 2 subgroups |
| info | `G01-S010-P0001-R002` | 3HYV | Listed in 2 subgroups |
| info | `G01-S010-P0003` | 2GMH | Listed in 2 subgroups (noted in description) |
| info | `G01-S010-P0003-R001` | 2GMJ | Listed in 2 subgroups |
| info | `G01-S011-P0001` | 2QCU | Listed in 2 subgroups (noted in description) |
| info | `G01-S011-P0001-R001` | 2R46 | Listed in 2 subgroups |
| info | `G01-S011-P0001-R002` | 2R45 | Listed in 2 subgroups |
| info | `G01-S011-P0001-R003` | 2R4J | Listed in 2 subgroups |
| info | `G01-S011-P0001-R004` | 2R4E | Listed in 2 subgroups |
| info | `G01-S016-P0003` | 3BF0 | Listed in 2 subgroups (noted in description) |
| info | `G01-S016-P0003-R001` | 3BEZ | Listed in 2 subgroups |
| info | `G01-S016-P0004` | 3RST | Listed in 2 subgroups (noted in description) |
| info | `G01-S016-P0004-D001` | 4KWB | Listed in 2 subgroups (noted in description) |
| info | `G02-S008-P0008` | 8W5J | Listed in 2 subgroups |
| info | `G03-S002-P0003-R004` | 1ZLL | Listed in 2 subgroups |
| info | `G03-S048-P0029` | 8JM9 | Listed in 2 subgroups |
| info | `G03-S048-P0029-R001` | 8JMA | Listed in 2 subgroups |
| info | `G03-S048-P0029-R002` | 8X82 | Listed in 2 subgroups |
| info | `G03-S048-P0029-R003` | 8X83 | Listed in 2 subgroups |
| info | `G03-S048-P0029-R004` | 8X84 | Listed in 2 subgroups |

<a id="duplicate_subgroup_name"></a>
## Subgroup names used more than once

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| info | `G01-S005` |  | Subgroup name 'Hydrolases' is used by 2 subgroups; identify subgroups by key |
| info | `G03-S136` |  | Subgroup name 'Hydrolases' is used by 2 subgroups; identify subgroups by key |

<a id="empty_subgroup"></a>
## Subgroups with no entries

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| info | `G03-S050` |  | Subgroup 'G PROTEIN-COUPLED RECEPTORS (GPCRs)' has no entries (header subgroup?) |

<a id="empty_tag"></a>
## Empty formatting tags

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| info | `G03-S133-P0002-D001` | 7MPA | bibliography[0].title: empty formatting tag |

<a id="invisible_characters"></a>
## Zero-width characters

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| info | `G03-S037-P0010-D001` | 8QZ1 | description: contains zero-width characters (removed for display) |
| info | `G03-S157-P0002-D001` | 8SOQ | description: contains zero-width characters (removed for display) |

<a id="pdb_unreleased"></a>
## PDB IDs not yet released

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| info | `G03-S051-P0080-D002-R003` | 9PX2 | Not yet released by the PDB |
| info | `G03-S051-P0080-D002-R004` | 9PX3 | Not yet released by the PDB |
| info | `G03-S051-P0080-D002-R005` | 9PX4 | Not yet released by the PDB |
| info | `G03-S051-P0095-D004` | 9PJK | Not yet released by the PDB |
| info | `G03-S102-P0015` | 7UUV | Not yet released by the PDB |
| info | `G03-S145-P0005` | 7ROW | Not yet released by the PDB |

<a id="related_not_described"></a>
## Related entries not cited in the description

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| info | `G02-S003-P0026` | 3VZT | Related 3VZU not cited in the description |
| info | `G02-S007-P0016-D008` | 7RI7 | Related 7RJ5 not cited in the description |
| info | `G02-S010-P0001-D001` | 6ZAZ | Related 6ZLT not cited in the description |
| info | `G02-S012-P0001-D002` | 6XFJ | Related 6XFL not cited in the description |
| info | `G03-S002-P0001` | 4P6J | Related 2MUZ not cited in the description |
| info | `G03-S002-P0003` | 6MQU | Related 1ZLL not cited in the description |
| info | `G03-S035-P0018` | 3AQP | Related 2RRN not cited in the description |
| info | `G03-S036-P0019-D002` | 9E62 | Related 9E63 not cited in the description |
| info | `G03-S037-P0001-D019` | 7MHR | Related 7MJT not cited in the description |
| info | `G03-S037-P0046` | 2WLJ | Related 1XL6, 1XL4 not cited in the description |
| info | `G03-S037-P0052-D002` | 7Y1L | Related 7Y1N not cited in the description |
| info | `G03-S041-P0048` | 8D1E | Related 8D1N not cited in the description |
| info | `G03-S041-P0073-D002` | 8TS2 | Related 8TS1 not cited in the description |
| info | `G03-S049-P0002` | 6SGW | Related 6SGZ not cited in the description |
| info | `G03-S051-P0038-D031` | 7IO0 | Related 7INB not cited in the description |
| info | `G03-S051-P0055-D002` | 7RKM | Related 7RKX not cited in the description |
| info | `G03-S051-P0133-D002` | 7RMG | Related 7RMI not cited in the description |
| info | `G03-S051-P0192-D003` | 9M88 | Related 9M8V not cited in the description |
| info | `G03-S114-P0044-D003` | 7SVD | Related 7SVR not cited in the description |
| info | `G03-S133-P0001-D020` | 5XA7 | Related 5XAA not cited in the description |
| info | `G03-S133-P0006-D006` | 7D91 | Related 7D94 not cited in the description |
| info | `G03-S133-P0022-D005` | 7VGI | Related 7VGJ not cited in the description |
| info | `G03-S139-P0001` | 2VPZ | Related 2VPY, 2VPW, 2VPX not cited in the description |
| info | `G03-S142-P0016` | 6ZKO | Related 6ZKA, 6ZKB not cited in the description |
| info | `G03-S151-P0013` | 7Y5E | Related 7Y7A not cited in the description |

<a id="resolution_multiple"></a>
## Resolution with several values

| Severity | Record key | PDB | Finding |
|---|---|---|---|
| info | `G03-S051-P0035` | 2R4R | resolution '3.4/3.7 Å' has several values; using 3.4 |
