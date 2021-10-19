# wes_ko_ukbb

Pipeline for assessing compound heterozygous variation in UKBB using HAIL and SAIGE.

## scripts

### 00_dbNSFP 
Use VEP to generate annotations for dbNSFP including REVEL_score and CADD_score. So far, we have been unable to change the change the hail.vep condig file to also use the dbSNP library. For this reason, we are generating these annotaitons externally, and them adding them later to hail matrix tables.

**input:**: quality controlled exome sequenced variants (.vcf.gz)
**output:**: VEP annotated variants using dbNSFP and LOFTEE (.vcf)

### 01_hail_vep
Use hl.vep and gnomad.process_consequence to create a hail table (of rows) that can be imported later onto phased or unphased data.

**input**: data/unphased/
**output**: data/vep/hail/ukb_200k_chr_vep.ht


