# Conditional analysis for nearby common variants

**Strategy:** "Is the the compound heterozyogus association signal driven by nearby common variants?" To ensure that this is not the case, we find nearby common variants that are significantly associated with our trait of interest. Then, we iteratively condition on all nearby variants in the locus until no variants are significantly associated with the trait (exome-wide). We append pseudo variants with this list of variants and re-run the primary analysis, but this time conditionion on the nearby common variants (and polygenic risk score). 


| Script                            | Description                                   | Requires output from |
| -                                 | -                                             | - |
| `00_get_gene_positions.R` | Get GRCh38 positions of genes on autosomes using R::biomaRt | N/A |
| `01_extract_intervals.sh` | 1. Load SAIGE results from primary analysis. </br> 2. Adjust SPA P-values with FDR. </br> Save markers (genes) with FDR<10%. | `00_get_gene_positions.R` |
| `02_filter_genotypes.sh` | 1. Load gene intervals and append 1MB upstream and downstream. </br> 2. Get UKB 500k genotype calls and filter to [QCed samples](https://github.com/astheeggeggs/SAIGE_gene_munging/tree/main/QC_scripts) </br> 3. Filter to variants with missing<10%, imputation INFO > 0.8 and with MAF>=1% or do phenotype specific MAF filtering. </br> 5. Filter to variants in gene intervals from step (1). </br> 6. Write out VCF. | `01_extract_intervals` |
| `03_spa_iter_common.sh`, `_spa_iter_common.sh` and `_spa_iter_common.R`| 1. Load filtered genoypes. </br> 2. Run conditional common analysis loop: <br> A. Run step 2 from SAIGE generating a list of condotional markers and resulting P-values. <br> B. Loop over markers and select most significant marker (based on oredering the T-statistics) <br> C. If it passes the SPA P-value threshold of P-value<=5e-6 (exome-wide), then add marker to conditional marker list. | `02_filter_genotypes.sh` |
| `04_extract_marker_gt.sh`, `04_extract_marker_gt.R` and `04_extract_marker_gt.py` | Based on the significant markers from the previous analysis, extract their corresponding genotypes across all chromosomes and save in a single VCF. | `03_spa_iter_common.sh` |
| `05_append_vcf_common.sh` and ` 05_append_vcf_common.py` | Append genotypes from common markers in previous step with pseudo markers from knockout analysis | `04_extract_marker_gt.sh` |
| `06_spa_cond_common.sh`, ` _spa_cond_common.sh` and `_spa_cond_common.R` | Run SAIGE on pseudo-variants while conditionion on common hits from previous steps. | `05_append_vcf_common.sh` |




