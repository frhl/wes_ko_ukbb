# Whole Exome Sequencing haplotype estimation


## Overview
| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `01_calls_gen.sh` | 1. Extract genotype calls data </br> 2. Perform liftover on variants. Subset to non-singletons variants, i.e. `MAC>=2`. Filter out variants variants that have more than 5% missingness.  | N/A |
| `01_wes_union_calls_gen.sh` | 1. Extract genotype calls data </br> 2. Perform liftover on variants. </br> 3. Extract quality controlled Whole Exome sequences (WES). </br> 4. Subset to samples present in both genotype calls and WES. </br> 5. Exclude any variants in calls, that are found among WES variants. </br> 6. Merge/combine the two datasets. </br> 7. Subset to non-singletons variants, i.e. `MAC>=2`. Filter out variants variants that have more than 5% missingness.  | N/A |

## Ressources

* [Shapeit4](https://odelaneau.github.io/shapeit4/)
* [Exome imputation pipeline](https://data.broadinstitute.org/lohlab/UKB_exomeWAS/code/imputation/) from Po-Ruh Loh


## references

* [Delaneau, O., Zagury, JF., Robinson, M.R. et al. Accurate, scalable and integrative haplotype estimation. Nat Commun 10, 5436 (2019)](https://www.nature.com/articles/s41467-019-13225-y)
* [Whole-exome imputation within UK Biobank powers rare coding variant association and fine-mapping analyses, Alison R. Barton, Maxwell A. Sherman, Ronen E. Mukamel, Po-Ru Loh](https://www.medrxiv.org/content/10.1101/2020.08.28.20180414v1)


