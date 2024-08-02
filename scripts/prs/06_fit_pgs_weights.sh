#!/usr/bin/env bash
#
# Fit weights for PGS
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=fit_pgs_weights
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/fit_pgs_weights_test.log
#SBATCH --error=logs/fit_pgs_weights_test.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task 10
#SBATCH --array=1-50
#SBATCH --qos="960_cpu"
# --begin=now+3hour
#

# Note: this script scales massively with parallelization
# 1 core takes 12-32h
# 6 cores tale <2h

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/06_fit_pgs_weights.R"
readonly rscript_check_prs="scripts/_check_prs_ok.R"

#readonly ldsc_dir="data/prs/ldsc"
readonly ldsc_dir="data/prs/ldsc_test"
readonly ld_dir="data/prs/hapmap/ld/matrix_unrel_kin"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/weights/auto_mle_false"

readonly ldsc_pvalue_cutoff="0.05"
readonly ldsc_n_eff_cutoff=100 # 20000
readonly ldpred_method="auto"

readonly cluster=$( get_current_cluster )
readonly index=$( get_array_task_id )

readonly file_binary="${pheno_dir}/dec22_phenotypes_binary_500k.tsv.gz"
readonly pheno_list_binary="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )

readonly file_cts="${pheno_dir}/curated_covar_phenotypes_cts_int_500k.txt.gz"
readonly pheno_list_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotype_cts=$( sed "${index}q;d" ${pheno_list_cts} )


export OPENBLAS_NUM_THREADS=1 # avoid two levels of parallelization

mkdir -p ${out_dir}

submit_ldpred() {
  set_up_ldpred2
  local phenotype=${1}
  local method=${2}
  local ldsc="${ldsc_dir}/ldsc_${phenotype}.rds"
  local out_prefix="${out_dir}/weights_${phenotype}"
  if [ ! -z ${phenotype} ]; then
     #if [ ! -f "${out_prefix}.txt.gz" ]; then
     if [ ! -f "${out_prefix}.rda" ]; then
        Rscript "${rscript}" \
          --ldsc "${ldsc}" \
          --ld_dir "${ld_dir}" \
          --method "${method}" \
          --ldsc_pvalue_cutoff "${ldsc_pvalue_cutoff}" \
          --ldsc_n_eff_cutoff "${ldsc_n_eff_cutoff}" \
          --out_prefix "${out_prefix}"
     else
        >&2 echo "${out_prefix}.rda already exists. Skipping.."
     fi
  fi
}

#submit_ldpred ${phenotype_binary} ${ldpred_method}
submit_ldpred ${phenotype_cts} ${ldpred_method}


