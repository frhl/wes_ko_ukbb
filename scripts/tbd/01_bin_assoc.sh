#!/bin/bash

# 22-02-15: run association for binary traits

#$ -cwd
#$ -N bin_assoc_logistic
#$ -o bin_assoc_logistic.log
#$ -e bin_assoc_logistic.errors.log
#$ -q short.qf
#$ -P lindgren.prjc
#$ -pe shmem 2

# in files
IN_SAMPLES="assoc_in_samples.tmp"
IN_GENES="assoc_in_genes.tmp"

module load R

## read in genes
cat ${IN_GENES}| while read GENE
do
   echo $GENE
   Rscript run_assoc_analysis.R "${PHENO}" "${IN_SAMPLES}" "${GENE}" "${SGE_TASK_ID}" "${OUTDIR}" "${COVARS}" "${MODEL}" "${KO_ENCODING}"
done

