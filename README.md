# wes_ko_ukbb


## 00_dbNSFP.sh 
Use VEP to generate annotations for dbNSFP including REVEL_score and CADD_score. So far, we have been unable to change the change the hail.vep condig file to also use the dbSNP library. For this reason, we are generating these annotaitons externally, and them adding them later to hail matrix tables.

## 01_hail_format
Run hail.vep and gnomad.process_consequence to generate two matrix tables: one for unphased singletons and another for phased non-singletons.  




