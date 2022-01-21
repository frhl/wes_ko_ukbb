#!/usr/bin/env bash
#
#$ -N summary_statistics
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/summary_statistics.log
#$ -e logs/summary_statistics.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1
#$ -V

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/sumstat"

readonly hail_script="scripts/prs/01_summary_statistics.py"

readonly covar_file="${pheno_dir}/covars1.csv"
readonly covariates=$( cat ${covar_file} )

readonly pheno_file="${pheno_dir}/curated_phenotypes.tsv" 
readonly pheno_list="${pheno_dir}/curated_phenotypes_cts_header.tsv"
readonly phenotype=$( cut -f${SGE_TASK_ID} ${pheno_list} )

readonly out_prefix="${out_dir}/test_ukb_imp_500k_${phenotype}"

#readonly chr=$( seq 1 22 | tr '\n' ' ' )
readonly chr="21"

mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy
module load OpenBLAS/0.3.8-GCC-9.2.0
export LD_PRELOAD=/apps/eb/skylake/software/OpenBLAS/0.3.8-GCC-9.2.0/lib/libopenblas.so
set -x
python3 "${hail_script}" \
   --chrom "${chr}" \
   --dataset "calls" \
   --phenotypes ${pheno_file} \
   --response ${phenotype} \
   --covariates ${covariates} \
   --liftover \
   --min_info 0.80 \
   --out_prefix "${out_prefix}"
set +x


