# Polygenic Risk Scoring using LDpred2
| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `00_bed_gen.sh` | 1. Load imputed genotypes or regular non-imputed genotypes. </br> 2. Subset to white brish europeans (based on Duncan's) QC. </br> 3.  </br> 3. Liftover to GRCh38 if needed. </br> 4. Write plink (.bed) output files.  | N/A |
| `01_summary_statistics.sh`| 1. Load imputed genotype or non-imputed genotypes. </br> 2. Subset to defined samples. </br> 3. Liftover to GRCh38 if needed. </br> 4. Load binary or continious phenotypes and covariates. </br> 5. Perform GWAS (Linear/logistic regression) on phenotype parallized by chromosome..|

