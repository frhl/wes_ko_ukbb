#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=fit_simulated_null
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/fit_simulated_null.log
#SBATCH --error=logs/fit_simulated_null.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --array=22
#SBATCH --begin=now+2hour


set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly array_idx=$( get_array_task_id )
readonly chr=$(get_chr ${array_idx} )

readonly plink_dir="data/saige/grm/input"
readonly grm_dir="data/saige/grm/input/dnanexus"
readonly covar_dir="data/phenotypes"
readonly pheno_dir="data/simulation/phenotypes_new"
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
    simulate_phenotypes ${1} ${2} ${3} ${4} ${5}
}

# main script
simulate_phenotypes() {

  local K=${1}
  local h2=${2}
  local b=${3}
  local pi=${4}
  local seed=${5}

  local prefix="ukb_wes_union_calls_h2_${h2}_b_${b}_pi_${pi}_K${K}_seed${seed}_chr${chr}"
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
   local jname="_fit_simulated_null"
   local lname="logs/_fit_simulated_null"
   local sim_job=$( sbatch \
          --account="${project}" \
          --job-name="${jname}" \
          --output="${lname}.log" \
          --error="${lname}.errors.log" \
          --chdir="$(pwd)" \
          --partition="${queue}" \
          --cpus-per-task="${nslots}" \
          --array="${tasks}" \
          --parsable \
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
           ${chr} )
 }

readonly project="lindgren.prj"
readonly queue="short"
readonly nslots="3"
readonly tasks="1-50"

#run_with_params 0.001 0.00 0.0 0.25 501
#run_with_params 0.005 0.00 0.0 0.25 502
#run_with_params 0.010 0.00 0.0 0.25 503
#run_with_params 0.020 0.00 0.0 0.25 504
#run_with_params 0.030 0.00 0.0 0.25 505
#run_with_params 0.050 0.00 0.0 0.25 506

#run_with_params 0.001 0.00 0.0 0.25 701
#run_with_params 0.003 0.00 0.0 0.25 702
#run_with_params 0.005 0.00 0.0 0.25 703
#run_with_params 0.008 0.00 0.0 0.25 704
#run_with_params 0.010 0.00 0.0 0.25 705
#run_with_params 0.020 0.00 0.0 0.25 706
#run_with_params 0.030 0.00 0.0 0.25 707
#run_with_params 0.050 0.00 0.0 0.25 708
#run_with_params 0.100 0.00 0.0 0.25 709

#run_with_params 0.001 0.00 0.0 0.25 601
#run_with_params 0.003 0.00 0.0 0.25 602
#run_with_params 0.005 0.00 0.0 0.25 603
#run_with_params 0.008 0.00 0.0 0.25 604
#run_with_params 0.010 0.00 0.0 0.25 605
#run_with_params 0.020 0.00 0.0 0.25 606
#run_with_params 0.030 0.00 0.0 0.25 607
#run_with_params 0.050 0.00 0.0 0.25 608
#run_with_params 0.100 0.00 0.0 0.25 609

run_with_params 0.001 0.00 0.5 0.25 401
run_with_params 0.001 0.01 0.5 0.25 401
run_with_params 0.001 0.02 0.5 0.25 401
run_with_params 0.001 0.05 0.5 0.25 401
run_with_params 0.001 0.10 0.5 0.25 401

run_with_params 0.001 0.00 1.0 0.25 401
run_with_params 0.001 0.01 1.0 0.25 401
run_with_params 0.001 0.02 1.0 0.25 401
run_with_params 0.001 0.05 1.0 0.25 401
run_with_params 0.001 0.10 1.0 0.25 401

run_with_params 0.001 0.00 2.5 0.25 401
run_with_params 0.001 0.01 2.5 0.25 401
run_with_params 0.001 0.02 2.5 0.25 401
run_with_params 0.001 0.05 2.5 0.25 401
run_with_params 0.001 0.10 2.5 0.25 401

run_with_params 0.001 0.00 10.0 0.25 401
run_with_params 0.001 0.01 10.0 0.25 401
run_with_params 0.001 0.02 10.0 0.25 401
run_with_params 0.001 0.05 10.0 0.25 401
run_with_params 0.001 0.10 10.0 0.25 401

