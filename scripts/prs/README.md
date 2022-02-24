# Polygenic Risk Scoring using LDpred2

Todo:
* Make holdout set
* Dan Howrigan used was 25/(2 x min(cases, controls)) 

## Overview
| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `00_bed_gen.sh` | 1. Load imputed genotypes or regular non-imputed genotypes. </br> 2. Subset to white-british europeans (Based on Duncan's QC) and Subset to unrelated individuals based on kinship coeffecient (Removing 1st, 2nd and 3rd degree relatives). </br> 3. Subset to HapMap3 SNPS with MAF >= 1% and filter out any genotypes with missigness > 5%. Also annotate with dbSNP. Liftover to GRCh38 </br> 4. Write plink (.bed)  output files.  | N/A |
| `01_sample_bed.sh` | 1. Load imputed genotypes or regular non-imputed genotypes. </br> 2. Subset to white-british europeans (Based on Duncan's QC) and subset to unrelated individuals. Randomly sample 10.000 individuals. </br> 3. Subset to HapMap3 SNPS with MAF >= 1% and filter out any genotypes with missigness > 5%. Also annotate with dbSNP. Liftover to GRCh38. </br> 4. Write plink (.bed)  output files.  | N/A |
| `02_merge_bed.sh` | 1. load plink (.bed) for all chromosomes available. </br> 2. Merge all chromosomes and write new combined plink (.bed) file.  | `01_sample_bed.sh` |
| `03_gwas.sh`| 1. Load unrelated white British samples. </br> 2. Load binary or continious phenotypes and covariates. Note: We only include binary phenotypes with at least 1250 cases. </br> 3. Do linear or logistic (wald) regression on the phenotypes to generate summary statistics. </br> 4. Submits merge script that will aggregate resulting GWAS into a single file. </br> 5. Write out summary statistics. | `00_bed_gen.sh` |
| `04_calc_ld.sh`| 1. Load plink (bed) containing 10.000 unrelated samples.  </br> 2. Calculate SNP-wise correlations for each chromosome. </br> 3. Save resulting correlation in an RDS-file.   | `01_sample_bed.sh` |
| `05_ldsc.sh`| 1. Load plink (bed) containing 10.000 unrelated samples. Also load GWAS and RDS file containing SNP-wise correlations. </br> 2. Perform [quality control](https://github.com/privefl/paper-ldpred2/blob/master/code/prepare-sumstats.R) on SNPs from GWAS. </br> 3. For each chromosome, open RDS file containg SNP-wise correlations. Create a sparse matrix file on the fly and calculate Linkage Disequilibrium (LD) for SNPs in GWAS that passed QC. </br> 4. Ensure that SNPs in GWAS, LD-matrix and correlation matrix has been aligned. </br> 5. run LDSC and obtain (genome-wide) heritability estimate for trait.  </br> 6. Save RDS-file with heritability estimate, QCed GWAS SNPs, and mapping file to SNPs used for LD-correlations (equivalent to GWAS file).     | `01_sample_bed.sh`, `03_gwas.sh` and `04_calc_ld.sh` |
| `06_prs.sh` | 1. Open full set of samples in plink (.bed) file. </br> 2. open RDS-file containing genome-wide h2 estimate. Partition the estimate with h2 * (n_chr/n_total) and use this as the initial estimate for heritability downstream. </br> 3. Match input genotypes with SNP-order of SNP-correlation matrix and GWAS SNPs. </br> 4. Open SNP-correlation matrix for current chromosome. </br> 4. Obtain polygenic risk scores with *LDPred2* using either the [infinitesimal](https://privefl.github.io/bigsnpr/reference/LDpred2.html) or [auto](https://privefl.github.io/bigsnpr/reference/LDpred2.html) mode.          | `05_ldsc.sh` and `00_bed_gen.sh`  |
| `07_aggr_prs.sh`| 1. Aggregate all PRS into a single file keyed by sample id.  | `06_prs.sh` |


## Ressources

* [LDpred2 tutorial](https://privefl.github.io/bigsnpr/articles/LDpred2.html) by Florian Privé
* [Anoter PRS tutorial](https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html) by Florian
* [Details and considerations](http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas) when conducting gwas in UK Biobank


## references

* [Florian Privé, Julyan Arbel, Bjarni J Vilhjálmsson, LDpred2: better, faster, stronger, Bioinformatics, Volume 36, Issue 22-2 2020](https://academic.oup.com/bioinformatics/article/36/22-23/5424/6039173)
* [Vilhjálmsson BJ, Yang J, Finucane HK, et al. Modeling Linkage Disequilibrium Increases Accuracy of Polygenic Risk Scores. Am J Hum Genet 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4596916/)


