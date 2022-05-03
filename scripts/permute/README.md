# Permutation pipeline

* Note, that `03_permute.sh` will result in the submission of thousands of jobs. Please test on chromosome 21 before submitting for all chromosomes.
* Note, that we currently consider phased and unphased hets as the same when specifying the underlying distibution for drawing knockouts (i.e. MAC is not relevant here).

## Overview
| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `01_calc_min_permutations.sh` | Calculate the minimum amount of permutations required based on P-values from primary knockout analysis  | `<primary spa analysis>`  |
| `02_count_hets.sh` | Count amount of phased, unphased hets and homozgygotes in knockout file. Note, this script assumes that the input has already been subsetted/filtered to damaging variants. The sampling distribution/probability for being a true knockout (pTKO) is determined by probability of all variants being on the same phase, i.e. pTKO = 1 - p(all on phase 1) - p(all on phase 2) = 1 - 2p(all on phase 1) `pTKO(k) = 1 - 2(1/2)^k` in whick k is th e number of heteterozygotes with `MAC>=1` in a gene.  | `<variant VCFs>`  |
| `03_permute.sh` | Master script stratified by chromosome. Will start a new script (`_permute.sh`) for each gene on the current chromosome. This script will then perform the following operations: </br>1. Check if any permutations are required based, and if so it will submit an array job (`_gene_permute`) with number of tasks equal to number of permutations currently required. If not, the script will exit and no more iterations will be performed. </br>2. Submit a array job (`_gene_spa.sh`) performing saddlepoint approximation with SAIGE for each phenotype that still requires permutations. </br>3. Submit a job (`_merge_spa.sh`) that merges the saddlepoint approximated phenotypes (one merged file for each phenotype). </br>4. Submit a job (`_calc_p.sh`) that calculates the current empirical P-value by comparing the true P-values and the shuffled P-values originating from the permuted phase. </br>5. Submit itself again but (`_permute.sh`) increasing the number of shuffles required by a factor fo 10 for phenotypes that are not yet appropiately permuted. | `02_count_hets.sh` and ``01_calc_min_permutations.sh``  |
| `04_collect.sh` | Aggregates the permuted output in a manner that can easily be read by other software.  | `03_permute.sh`  |


