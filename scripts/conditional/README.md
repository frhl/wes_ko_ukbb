# Conditional analysis pipeline
| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `00_get_gene_positions.R` | Get GRCh38 positions of genes on autosomes using R::biomaRt | N/A |
| `01_extract_intervals.sh` | 1. Load SAIGE results from primary analysis. </br> 2. Adjust SPA P-values with FDR. </br> Save markers (genes) with FDR<10%. | `00_get_gene_positions.R` |
| `02_filter_genotypes.sh` | 1. Load gene intervals and append 1MB upstream downstream. </br> 2. Get UKB 500k genotype calls and filter to [QCed samples](https://github.com/astheeggeggs/SAIGE_gene_munging/tree/main/QC_scripts) </br> 3. Filter to variants with missing<10% </br> 4. Filter to variants with MAF>=1% </br> 5. Filter to variants that are within gene intervals. </br> 6. Write out VCF. | `01_extract_intervals` |
| `03_spa_conditional.sh` | 1. Load genotype VCFs and apply SPA using the null models from primary analysis. </br> 2. Loop over markers and select most significant marker if it passes the SPA P-value threshold of P-value<=10e-6. Condition on this marker and repeat analysis. Restart loop until no markers are significant at the threshold given. </br> 3.Save conditioning markers. | `02_filter_genotypes.sh` |
| `04_merge_markers.sh` | 1. Merge markers (genes) from primary analysis with markers (variants) from previous contioning analysis.  | `03_spa_conditional.sh` |




