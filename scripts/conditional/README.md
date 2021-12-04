# Conditional analysis pipeline
| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `00_get_gene_positions.R` | Get GRCh38 positions of genes on autosomes using R::biomaRt | N/A |
| `01_extract_intervals.sh` |  | `00_get_gene_positions.R` |
| `02_filter_genotypes.sh` | 1. Load gene intervals and append 1MB upstream downstream. </br> 2. Get UKB 500k genotype calls and filter to [QCed samples](https://github.com/astheeggeggs/SAIGE_gene_munging/tree/main/QC_scripts) </br> 4. Filter to variants with missing<10% </br> 4. Filter to variants with MAF>=1% </br> 5. Filter to variants that are within gene intervals. </br> 6. Write out VCF. | `01_extract_intervals` |
| `03_spa_conditional.sh` |  | N/A |




