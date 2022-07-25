# Conditional analysis pipeline

Hi everyone,
Here's where I think we got to - please correct me if I missed anything!
* compount_het variants = aggregated variant (presence/absence of compound het,  where each variant in the compound het has MAF < 0.05), one for each gene.
* common variants = 0,1,2 encoded common variants in the region.
* rare variants = 0,1,2 encoded rare variants in the gene
* LOCO PRS = off chromosome PRS estimated from common variants using whatever PRS software we fancy
1. Filter to genes with a significant compound_het variant association, move to step 2.
2. Filter to genes with significant compound_het variant association after also conditioning on LOCO PRS, move to step 3.
3. Filter to genes with significant compound_het variant association after additionally conditioning on all common variant signals in the region using an iterative procedure. Are there any left?
4. We now have a couple of options:
a. What Wei called 'reciprocal conditioning'. Run standard SAIGE gene on rare variants in these genes conditioning on off chr PRS, common signal nearby - get a set of P-values. Run the same test, but, additionally, condition on the compound het signal - did the signal disappear for any of the tests? If it did, we have a candidate association that's driven by compound hets.
b. Just run the brute force approach of conditioning out all rare variant PTVs in the gene, and see if there's still a significant association between presence/absence of the compound het.

## Overview

* `common`: Condition on nearby common variants
* `rare`: Condition on rare variants in gene
* `combined`: Either `reciprocal` with conditioning on common variants, prs, and knockouts or `brute_force`: Condition on common variants, prs and rare variants



