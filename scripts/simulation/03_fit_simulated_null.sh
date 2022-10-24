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
readonly grm_dir="data/saige/grm/input/dnanexus"
readonly covar_dir="data/phenotypes"
readonly pheno_dir="data/simulation/phenotypes"
readonly out_dir="data/simulation/saige/step1/binary"

readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly plink_file="${grm_dir}/ukb_eur_200k_grm_grch38_rv_merged"
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

  #local trait_type="quantitative"
  local trait_type="binary"
  #local inv_normalize="TRUE"
  local inv_normalize="FALSE"
  local phenotype="case" # or "case" for binary
  #local phenotype="y"
  #local phenotype="case"
  fit_null
}

fit_null() {
   local out_prefix="${out_dir}/${prefix}_${phenotype}"
   set -x
   qsub -N "_sim_fit_null${SGE_TASK_ID}" \
     -q "short.qe" \
     -pe shmem 2 \
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


# this works and results in nice phenotypes
#run_with_params 0.00 0.00 0.00 0.01 0.01 600
#run_with_params 0.001 0.10 99.0 0.20 0.20 601
#run_with_params 0.002 0.10 99.0 0.20 0.20 602
#run_with_params 0.005 0.10 99.0 0.20 0.20 603
#run_with_params 0.01 0.10 99.0 0.20 0.20 604
#run_with_params 0.02 0.10 99.0 0.20 0.20 605
#run_with_params 0.05 0.10 99.0 0.20 0.20 606

#run_with_params 0.001 0.10 99.0 1.00 1.00 601
#run_with_params 0.002 0.10 99.0 1.00 1.00 602
#run_with_params 0.005 0.10 99.0 1.00 1.00 603
#run_with_params 0.01 0.10 99.0 1.00 1.00 604
#run_with_params 0.02 0.10 99.0 1.00 1.00 605
#run_with_params 0.05 0.10 99.0 1.00 1.00 606


run_with_params 0.001 10.0 0.10 0.20 0.20 501
run_with_params 0.001 0.10 0.10 0.20 0.20 501
run_with_params 0.001 0.10 1.00 0.20 0.20 502
run_with_params 0.001 0.10 10.0 0.20 0.20 503
run_with_params 0.001 0.10 99.0 0.20 0.20 504

#run_with_params 0.005 10.0 0.10 0.20 0.20 501
#run_with_params 0.005 0.10 0.10 0.20 0.20 501
#run_with_params 0.005 0.10 1.00 0.20 0.20 502
#run_with_params 0.005 0.10 10.0 0.20 0.20 503
#run_with_params 0.005 0.10 99.0 0.20 0.20 504

#run_with_params 0.01 10.0 0.10 0.20 0.20 501
#run_with_params 0.01 0.10 0.10 0.20 0.20 501
#run_with_params 0.01 0.10 1.00 0.20 0.20 502
#run_with_params 0.01 0.10 10.0 0.20 0.20 503
#run_with_params 0.01 0.10 99.0 0.20 0.20 504

#run_with_params 0.10 10.0 0.10 0.20 0.20 501
#run_with_params 0.10 0.10 0.10 0.20 0.20 501
#run_with_params 0.10 0.10 1.00 0.20 0.20 502
#run_with_params 0.10 0.10 10.0 0.20 0.20 503
#run_with_params 0.10 0.10 99.0 0.20 0.20 504

