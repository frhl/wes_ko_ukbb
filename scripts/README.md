# Compound Heterozygous pLOF pipeline


## Overview

* `scripts`: Main scripts
* `scripts/conditional`: Identify common and rare variants for conditional analysis
* `scripts/permute`: Permute genetic phase to empirical P-values across phentoypes
* `scripts/phasing`: Phasing whole exome sequencing data
* `scripts/post_hoc`: Summary of knockouts and variants
* `scripts/prs`: Polygenic risk scores using ldpred2
* `scripts/saige_gene`: set-based analysis with SAIGE-GENE+
* `scripts/simulation`: Simulation of compound heterozygous effects
* `scripts/survival`: Generate tables for surivival analysis


## Main scripts
| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `00_phenotypes.sh` | Filter phenotype table based on specific samples and sex. Count up cases and controls for binary phenotypes. Export header of phenotypes.  | N/A |
| `01_hail_vep.sh` | 1. Generate VEP annotations with gnomAD's `process_consequence` </br> 2. Export as hail tables. | N/A |
| `02_union_mts.sh` | Combine unphased variants with `AC==1` (singletons) with fully phased variants (non-singletons) from the phasing pipeline. Resulting MatrixTable contains phased and unphased gentoypes.   | [phasing pipeline](https://github.com/frhl/wes_ko_ukbb/tree/critter_tunes/scripts/phasing) |
| `03_annotate.sh` | 1. Perform subsets based Duncan's QCed samples and variants. </br> 2. Remove samples that have withdrawn consent. </br> 3. Annotate with VEP. </br> 4. Write MatrixTable of phased data.     | `02_union_mts.sh`</br>`01_hail_vep.sh` |
| `04_export_csqs.sh` | Export VEP consequences in a manner compatible with SAIGE-GENE+ (argument `groupFile`) to be used downstream. Also used for post-hoc interrogation of knockout variants. | `01_hail_vep.sh`</br>`03_annotate.sh` |
| `05_knockouts.sh` | 1. Read in phased (and optionally unphased) genotypes. </br> 2. Collapse variants by [consequence annotations](https://github.com/frhl/wes_ko_ukbb/blob/main/utils/modules/python/ko_utils/analysis.py) and select most deleterious candidate variant by gene. </br> 3. For each consequence category, filter to variants in the category. </br> 4. Using both unphased singletons and phased data, combine dosage based on </br> &emsp; i. There's n phased PTVs already in the gene, with at least 2 in trans, so we're done (it's already a knockout). </br> &emsp; ii. There's n phased PTVs already in the gene, all in cis, and k unphased singleton PTVs. We assume equal probability on being on each haplotype, therefore P(knockout) = 1 - (1/2)^k. </br> &emsp; iii. There's 0 phased PTVs already in the gene, and k >= 2 unphased singleton PTVs. In which case P(knockout) = (1 - 2*(1/2)^k). </br> 5. Write a synthethic VCF where rows correspond to gene, columns correspond to sample, and entries correspond to dosage DT), i.e., probability of being a knockout multiplied by 2.   | `03_annotate.sh` |
| `06_grm.sh` | 1. Load all autosomes from ukbiobank genotype calls and subset to variants in `/well/lindgren/UKBIOBANK/DATA/QC/ukb_snp_qc.txt` with `in_Relatedness==1`.  </br> 2. Perform sample subsets based Duncan's [QCed samples](https://github.com/astheeggeggs/SAIGE_gene_munging/tree/main/QC_scripts). </br> 3. Write plink (.bed, .bim, .fam) files. </br>  4. use SAIGE to fit sparse genetic relatedness matrix.  | N/A |
| `07_spa_null.sh` | Use SAIGE to fit a null model for each phenotype. | `00_phenotypes.sh` </br> `05_create_grm.sh` |
| `08_spa_test.sh` | Uuse SAIGE to fit a null model for each phenotype. Note that `min_mac` should be set. Wei recommends at least 4. </br> 2. For each trait `_spa_test.sh` will be called and resulting SPA P-values and test staistic returned.  | `00_phenotypes.sh`</br>`04_knockouts.sh`</br>`05_create_grm.sh`</br>`06_spa_null.sh` |





