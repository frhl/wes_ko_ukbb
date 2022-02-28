# Whole Exome Sequencing haplotype estimation

## Overview
| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `01_prefilter_wes.sh` | 1. Load QCed WES data </br> 2. Subset to non-singletons variants, i.e. `MAC>=2`. Filter out variants variants that have more than 5% missingness.  | N/A |
| `02_calls_gen.sh` | 1. Extract genotype calls data </br> 2. Perform liftover on variants. Subset to non-singletons variants, i.e. `MAC>=2`. Filter out variants variants that have more than 5% missingness.  | N/A |
| `02_wes_union_calls_gen.sh` | 1. Extract genotype calls data </br> 2. Perform liftover on variants. </br> 3. Extract quality controlled Whole Exome sequences (WES). </br> 4. Subset to samples present in both genotype calls and WES. </br> 5. Exclude any variants in calls, that are found among WES variants. </br> 6. Merge/combine the two datasets. </br> 7. Subset to non-singletons variants, i.e. `MAC>=2`. Filter out variants variants that have more than 5% missingness.  | N/A |
| `03_wes_naive_phasing.sh` | 1. Open whole exome sequencing data </br> 2. Run SHAPEIT4 phasing with [default parameters](https://odelaneau.github.io/shapeit4/) using the `--sequencing` flag. </br> 3. Index result VCF and create `.tbi` files. </br> 4. Calculate combined (trio) switch errors using BCFtools. | `01_prefilter_wes.sh` |
| `03_calls_naive_phasing.sh` | 1. Open genotype calls </br> 2. Run SHAPEIT4 phasing using default parameters. </br> 3. Index result VCF and create `.tbi` files. </br> 4. Calculate combined (trio) switch errors using BCFtools. | `02_calls_gen.sh` |
| `03_wes_union_calls_phasing.sh` | 1. Open WES+CALLS data. </br> 2. Run SHAPEIT4 phasing with [default parameters](https://odelaneau.github.io/shapeit4/) using the `--sequencing` flag. </br> 3. Index result VCF and create `.tbi` files. </br> 4. Calculate combined (trio) switch errors using BCFtools. | `01_wes_union_calls_gen.sh` |
| `03_calls_naive_ser.sh`, `03_wes_naive_ser.sh` and `03_wes_union_calls_ser.sh` | 1. Calculate trio switch error rates on indidiual sites using the [modified](https://github.com/frhl/wes_ko_ukbb/tree/main/utils/bcftools/plugins) BCFtools plugin. Calculate variant level statistics using the `bcftools +fill-tags` command. </br> 2. Combine variant-level statistics and switch errors for each site. Write out a switch error summary file stratified by variant position and individual. | `02_wes_union_calls_phasing.sh` |
| `04_phase_chunks.sh` | Using either `SHAPEIT4 or `Eagle2` to phase data chromosome-wise in chunks ensuring overlap between the sets. The sizes are controlled by the `phasing_region_overlap`, `phasing_region_overlap` and `max_phasing_region_size` parameters.   | `02_wes_union_calls_gen.sh` |
| `05_merge_chunks.sh` | Subset each phased chunk to the flank that correspond to half the length of the overlap. | `04_phase_chunks.sh` |


## Ressources

* [SHAPEIT4](https://odelaneau.github.io/shapeit4/) manual
* [Exome imputation pipeline](https://data.broadinstitute.org/lohlab/UKB_exomeWAS/code/imputation/) from Po-Ru Loh
* [BCFtools](https://samtools.github.io/bcftools/) version 1.12

## References

* [Delaneau, O., Zagury, JF., Robinson, M.R. et al. Accurate, scalable and integrative haplotype estimation. Nat Commun 10, 5436 (2019)](https://www.nature.com/articles/s41467-019-13225-y)
* [Whole-exome imputation within UK Biobank powers rare coding variant association and fine-mapping analyses, Alison R. Barton, Maxwell A. Sherman, Ronen E. Mukamel, Po-Ru Loh](https://www.medrxiv.org/content/10.1101/2020.08.28.20180414v1)


