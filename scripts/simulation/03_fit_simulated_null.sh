#!/usr/bin/env bash
#
#$ -N fit_simulated_null
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/fit_simulated_null.log
#$ -e logs/fit_simulated_null.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 22
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly chr="${SGE_TASK_ID}"
readonly plink_dir="data/saige/grm/input"
readonly grm_dir="data/saige/grm/input"
readonly covar_dir="data/phenotypes"
readonly pheno_dir="data/simulation/phenotypes"
readonly out_dir="data/simulation/saige/step1/cts"

readonly grm_mtx="${grm_dir}/211102_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly plink_file="${plink_dir}/chunks/ukb_wes_200k_sparse_autosomes_mrg"
readonly covar_file="${covar_dir}/covars2.csv"
readonly covariates=$( cat ${covar_file} )

readonly bash_script="scripts/simulation/_fit_simulated_null.sh"
readonly rscript="utils/saige/step1_fitNULLGLMM.R"

mkdir -p ${out_dir}


# wrapper for running main script
run_with_params() {
    simulate_phenotypes ${1} ${2} ${3} ${4} ${5} ${6}
}

# main script
simulate_phenotypes() {

  local K=0.1
  local h2=${1}
  local var_beta=${2}
  local var_theta=${3}
  local pi_beta=${4}
  local pi_theta=${5}
  local seed=${6}

  local vars="${var_beta}_${var_theta}"
  local pis="${pi_beta}_${pi_theta}"

  local prefix="ukb_eur_h2_${h2}_var_${vars}_pi_${pis}_K${K}_seed${seed}_chr${chr}"
  #local prefix="ukb_eur_h2_${h2s}_pi_${pis}_K${K}_seed${seed}_chr${chr}"
  local pheno_file="${pheno_dir}/${prefix}_phenos.tsv.gz"

  local trait_type="quantitative"
  #local trait_type="binary"
  local inv_normalize="TRUE"
  #local inv_normalize="FALSE"
  #local phenotype="y" # or "case" for binary
  local phenotype="y"
  #local phenotype="case"
  fit_null
}

fit_null() {
   local out_prefix="${out_dir}/${prefix}_${phenotype}"
   set -x
   qsub -N "_sim_fit_null${SGE_TASK_ID}" \
     -q "short.qe" \
     -pe shmem 1 \
     -t ${tasks} \
     ${bash_script} \
     ${rscript} \
     ${plink_file} \
     ${pheno_file} \
     ${phenotype} \
     ${covariates} \
     ${trait_type} \
     ${inv_normalize} \
     ${out_prefix} \
     ${grm_mtx} \
     ${grm_sam} \
     ${chr}
   set +x
 }

readonly tasks="1-5"

run_with_params 0.002 0.10 0.10 0.20 0.20 201
run_with_params 0.002 5.00 0.10 0.20 0.20 201
run_with_params 0.002 0.10 5.00 0.20 0.20 201

run_with_params 0.003 0.10 0.10 0.20 0.20 201
run_with_params 0.003 5.00 0.10 0.20 0.20 201
run_with_params 0.003 0.10 5.00 0.20 0.20 201

run_with_params 0.005 0.10 0.10 0.20 0.20 201
run_with_params 0.005 5.00 0.10 0.20 0.20 201
run_with_params 0.005 0.10 5.00 0.20 0.20 201



