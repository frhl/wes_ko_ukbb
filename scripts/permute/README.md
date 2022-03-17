# Permutation pipeline

## Overview
| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `01_calc_min_permutations.sh` | Calculate the minimum amount of permutations required based on P-values from primary knockout analysis  | `<primary spa analysis>`  |
| `02_count_hets.sh` | Count amount of phased, unphased hets and homozgygotes in knockout file (assuming subsetted to damaging variants). The sampling distribution (pTKO) for being a true knockout is based on the phased het counts.  | `<variant VCFs>`  |
| `03_array_permute.sh` | For each gene, realize actual knockouts based on probability of being a true knockout, and create probabilistic encoding (based on non-phased hets), and write to vcf. Replicate with respect to permutations needed.    | `02_count_hets.sh` and ``01_calc_min_permutations.sh``  |
| `04_array_spa.sh` | For each permuted gene (and its replicates), perform saddle point approximation to generate empirical null.    | `02_count_hets.sh`  |



