#!/bin/bash

#$ -cwd
# -q short.qf
#$ -P lindgren.prjc
# -pe shmem 2

# in files
#IN_SAMPLES="assoc_in_samples.tmp"
#IN_GENES="assoc_in_genes.tmp"

module load R

## read in genes
cat ${GENES}| while read GENE
do
   Rscript run_assoc_analysis.R "${PHENO}" "${SAMPLES}" "${GENE}" "${SGE_TASK_ID}" "${OUTDIR}" "${COVARS}" "${MODEL}" "${KO_ENCODING}"
done




