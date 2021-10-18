# wes_ko_ukbb

Pipeline for assessing compound heterozygous variation in UKBB using HAIL and SAIGE.

## scripts

### 00_dbNSFP.sh 
Use VEP to generate annotations for dbNSFP including REVEL_score and CADD_score. So far, we have been unable to change the change the hail.vep condig file to also use the dbSNP library. For this reason, we are generating these annotaitons externally, and them adding them later to hail matrix tables.

**input:** 
* quality controlled exome sequenced variants (.vcf.gz)

**output:**
* VEP annotated variants using dbNSFP and LOFTEE (.vcf)

### 01_hail_format.sh
Run hail.vep and gnomad.process_consequence to generate two matrix tables: one for unphased singletons and another for phased non-singletons. Results from dbNSFP will me merged into the matrix table. The table entries will also be re-annotated with DP and GQ, so that they can be used for filtering in downstream analysis.   

**input:**
* quality controlled exome sequenced variants (.vcf.gz or .mt/)
* phased data (.vcf.gz)
* VEP annotated variants (.vcf)

**output:**
* matrix table of unphased singletons
* matrix table of phased non-singletons

### 02_knockouts.sh
Performs various filtering operations on matrix tables based on variant and sample-level metrics. We consider three categories of mutations that are manually specified: PTV+damaging_missense, PTV and synonomous as a negative control. Then, it will aggregate variants by genes and determine what individuals that harbor homozygous, heterozygous or compound heterozygous variations. Additionally, the script will output a synthethic VCF of markers (genes) in which the dosage (DT) has been encoded as P(knockout).

**input:**
* matrix table of unphased singletons
* matrix table of phased non-singletons

**output:**
* burden matrix:
* knockout_matrix: 
* saige_vcf:


### 03_create_grm.sh


### 04


### 05 



