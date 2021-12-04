#!/usr/bin/env bash


#$ -cwd
#$ -N phase
#$ -o phase.log
#$ -e phase.errors.log
#$ -q long.qf 
#$ -pe shmem 20
#$ -V
#$ -P lindgren.prjc
#$ -t 1-22


## Variants and individuals
INDIVIDUALS="ALL" #80000
CHR=${SGE_TASK_ID}

## IN_FILES

MAP_FILE="/well/lindgren/flassen/software/SHAPEIT4/b37.gmap/chr${CHR}.b37.gmap.gz"
REF="/well/1000G/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"
#REF="/well/1000G/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr${CHR}.phase1_intergrated_calls.20101123.snps_indels_svs.genotypes.vcf.gz"


module load SHAPEIT4/4.1.3-foss-2019b
shapeit4 --input ${OUT_ROOT}${IN_VCF} \
  --map ${MAP_FILE} \
  --region ${CHR} \
  --thread 40 \
  --reference ${REF} \
  --output ${OUT_ROOT}${OUT_PHASED}



