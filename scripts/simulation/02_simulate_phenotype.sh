#!/usr/bin/env bash
#SBATCH --account=lindgren.prj
#SBATCH --job-name=simulate_phenotype
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/simulate_phenotype.log
#SBATCH --error=logs/simulate_phenotype.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --array=22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly bash_script="scripts/simulation/_simulate_phenotype.sh"
readonly merge_script="scripts/simulation/_merge_phenotype.sh"

readonly pheno_dir="data/phenotypes"
readonly pheno_file="${pheno_dir}/filtered_covar_phenotypes_cts.tsv.gz"

readonly covar_file="${pheno_dir}/covars2.csv"
readonly covariates=$( cat ${covar_file} )

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )
readonly in_dir="data/simulation/knockouts"
readonly in_prefix="${in_dir}/ukb_wes_union_calls_100k_encoded_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/simulation/phenotypes_new"

mkdir -p ${out_dir}

run_with_params() {
  simulate_phenotypes ${1} ${2} ${3} ${4} ${5}
}

simulate_phenotypes() {

  local K=${1}
  local h2=${2}
  local b=${3}
  local pi=${4}
  local seed=${5}

  local out_prefix="${out_dir}/ukb_wes_union_calls_h2_${h2}_b_${b}_pi_${pi}_K${K}_seed${seed}_chr${chr}"
  local out_phenotypes="${out_prefix}_phenos.tsv.gz"
  
  local jname="_simulate_phenotypes"
  local lname="logs/_simulate_phenotypes"
  local sim_job=$( sbatch \
        --account="${project}" \
        --job-name="${jname}" \
        --output="${lname}.log" \
        --error="${lname}.errors.log" \
        --chdir="$(pwd)" \
        --partition="${queue}" \
        --cpus-per-task="${nslots}" \
        --array=${tasks} \
        --open-mode="append" \
        --parsable \
        "${bash_script}" \
        "${in_prefix}" \
        "${in_type}" \
        "${h2}" \
        "${b}" \
        "${pi}" \
        "${K}" \
        "${seed}" \
        "${out_prefix}" )

  local jname="_merge_phenotypes"
  local lname="logs/_merge_phenotypes"
  local merge_job=$( sbatch \
        --account="${project}" \
        --job-name="${jname}" \
        --output="${lname}.log" \
        --error="${lname}.errors.log" \
        --chdir="$(pwd)" \
        --partition="${queue}" \
        --cpus-per-task="${nslots}" \
        --open-mode="append" \
        --parsable \
        --dependency=afterok:${sim_job} \
        "${merge_script}" \
        "${out_prefix}" \
        "${pheno_file}" \
        "${covariates}" \
        "${out_phenotypes}" )

}

###############
# main script #
###############


readonly project="lindgren.prj"
readonly queue="short"
readonly nslots="2"
readonly tasks=1-20

# K, h2, b, pi, seed
#run_with_params 0.001 0.00 0.0 0.25 701
#run_with_params 0.003 0.00 0.0 0.25 702
#run_with_params 0.005 0.00 0.0 0.25 703
#run_with_params 0.008 0.00 0.0 0.25 704
#run_with_params 0.010 0.00 0.0 0.25 705
#run_with_params 0.020 0.00 0.0 0.25 706
#run_with_params 0.030 0.00 0.0 0.25 707
#run_with_params 0.050 0.00 0.0 0.25 708
#run_with_params 0.100 0.00 0.0 0.25 709

#run_with_params 0.005 0.00 0.0 0.25 401
#run_with_params 0.005 0.01 0.0 0.25 401
#run_with_params 0.005 0.02 0.0 0.25 401
#run_with_params 0.005 0.05 0.0 0.25 401
#run_with_params 0.005 0.10 0.0 0.25 401

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








