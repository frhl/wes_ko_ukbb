# Compound Heterozygous pLOF pipeline

| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `00_dbNSFP.sh` | 1. Run VEP with dbNSFP plugin and generate REVEL_score and CADD_phred annotations. </br> 2. Exported as vcf.  | N/A |
| `01_hail_vep.sh` | 1. Generate VEP annotations with gnomAD's `process_consequence` and combine with previous dbNSFP annotations </br> 2. Export as hail tables. | `00_dbNSFP.sh` |
| `02_qc_input_mt.sh` | 1. Perform subsets based Duncan's [QCed samples and variants](https://github.com/astheeggeggs/SAIGE_gene_munging/tree/main/QC_scripts). </br> 2. Annotate with sites found in gnomAD. </br> 3. Annotate with sites found in imputed data. </br> 4. Export variant and sample level statistics using `hl.sample_qc` and `hl.variant_qc` as .tsv files. </br> 5. Annotate and export VEP consequence counts by gene. | `01_hail_vep.sh` |
| `03_create_mt.sh` | 1. Perform subsets based Duncan's QCed samples and variants. </br> 2. Annotate with VEP. </br> 3. Write MatrixTable of phased data. </br> 4. Subset unphased WES data with `AC==1` to get singletons and export as MatrixTable.    | `01_hail_vep.sh` |
| `04_knockouts.sh` | 1. Read in phased and unphased (singletons) simultanously. </br> 2. Collapse variant consequences into the [4-5 categories](https://github.com/frhl/wes_ko_ukbb/blob/main/utils/modules/python/ko_utils/analysis.py) and select most deleterious candidate variant by gene. </br> 3. For each consequence category, filter to variants in the category. </br> 4. Using both unphased singletons and phased data, combine dosage based on </br> &emsp; i. There's n phased PTVs already in the gene, with at least 2 in trans, so we're done (it's already a knockout). </br> &emsp; ii. There's n phased PTVs already in the gene, all in cis, and k unphased singleton PTVs. We assume equal probability on being on each haplotype, therefore P(knockout) = 1 - (1/2)^k. </br> &emsp; iii. There's 0 phased PTVs already in the gene, and k >= 2 unphased singleton PTVs. In which case P(knockout) = (1 - 2*(1/2)^k). </br> 5. Write a synthethic VCF where rows correspond to gene, columns correspond to sample, and entries correspond to dosage DT), i.e., probability of being a knockout multiplied by 2.   | `03_create_mt.sh` |
| `05_make_tabix.sh` | 1. Use BCFtools to make `.csi` indexing for all .vcf.bgz files.   | `04_knockouts.sh` |
| `06_create_grm.sh` | 1. Load all autosomes from ukbiobank imputed data autosomes and to variants in `/well/lindgren/UKBIOBANK/DATA/QC/ukb_snp_qc.txt` with `in_Relatedness==1`.  </br> 2. Perform subsets based Duncan's [QCed samples](https://github.com/astheeggeggs/SAIGE_gene_munging/tree/main/QC_scripts). </br> 3. Write plink files. </br>  4. use SAIGE to fit sparse genetic relatedness matrix  | N/A |
| `07_fit_null_saige_binary.sh` | 1. use SAIGE to fit a null model for each binary phenotype | `06_create_grm.sh` |
| `07_fit_null_saige_cts.sh` | 1. use SAIGE to fit a null model for each continuous phenotype | `06_create_grm.sh` |
| `08_spa_test.sh` | 1. use SAIGE to fit a null model for each continuous phenotype. </br> 2. For each trait (based on SGE_TASK_ID index), `_spa_test.sh` will be called.  | `07_fit_null_saige_binary.sh` </br> `04_knockouts.sh` </br> `05_make_tabix.sh` |

# Summary



