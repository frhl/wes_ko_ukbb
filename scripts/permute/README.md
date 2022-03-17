# Permutation pipeline

## Overview
| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `01_calc_min_permutations.sh` | Calculate the minimum amount of permutations required based on P-values from primary knockout analysis  | `<primary knockout analysis>`  |
| `02_count_hets.sh` | Count amount of phased, unphased hets and homozgygotes in knockout file (assuming subsetted to damaging variants). The probability of being true knockout is calculated based on these numbers (E.g. 2 phased hets, result in 50% chance of knockout).  | `NA`  |
| `03_array_permute.sh` | For each gene, realize actual knockouts based on probability of being a true knockout, and create probabilistic encoding (based on non-phased hets), and write to vcf. Repeat depending on permutations needed.    | `02_count_hets.sh` and ``01_calc_min_permutations.sh``  |
| `04_array_spa.sh` | For each permuted gene (and its replicates), perform saddle point approximation to generate empiricall null.    | `02_count_hets.sh`  |





