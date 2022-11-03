#!/usr/bin/env bash
#
#$ -N simulate_phenotype
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/simulate_phenotype.log
#$ -e logs/simulate_phenotype.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 22
#$ -V

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

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly in_dir="data/simulation/knockouts"
readonly in_prefix="${in_dir}/ukb_eur_100k_knockout_recessive_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/simulation/phenotypes"

mkdir -p ${out_dir}

run_with_params() {
  simulate_phenotypes ${1} ${2} ${3} ${4} ${5}
}

simulate_phenotypes() {

  local K=0.1
  local h2=${1}
  local var_beta=${2}
  local var_theta=${3}
  local pi=${4}
  local seed=${5}

  local vars="${var_beta}_${var_theta}"

  local out_prefix="${out_dir}/ukb_eur_h2_${h2}_var_${vars}_pi_${pi}_K${K}_seed${seed}_chr${chr}"
  local out_phenotypes="${out_prefix}_phenos.tsv.gz"
  
  local sim_name="sim_c${SGE_TASK_ID}_h2_${h2}_var_${vars}_pi_${pi}"
  local mrg_name="mrg_c${SGE_TASK_ID}_h2_${h2}_var_${vars}_pi_${pi}"

  set -x
  if [ 1 -eq 1 ]; then
  qsub -N "${sim_name}" \
       -t ${tasks} \
      -q ${queue} \
       -pe shmem ${nslots} \
       ${bash_script} \
       ${in_prefix}\
       ${in_type} \
       ${h2} \
       ${var_beta} \
       ${var_theta} \
       ${pi} \
       ${K} \
       ${seed} \
       ${out_prefix}
  fi

  qsub -N "${mrg_name}" \
       -hold_jid ${sim_name} \
       -t 1 \
       -q "short.qc" \
       -pe shmem 2 \
       ${merge_script} \
       ${out_prefix}\
       ${pheno_file} \
       ${covariates} \
       ${out_phenotypes}
  set +x
}

###############
# main script #
###############

readonly queue="short.qc"
readonly nslots="2"
readonly tasks=1-5 #-2

#simulate_phenotypes 0.00 0.00 0.00 0.01 0.01

#run_with_params 0.00 0.00 0.00 0.01 0.01 600
#run_with_params 0.01 95.00 0.01 0.25 0.25 601
#run_with_params 0.05 95.00 0.01 0.25 0.25 602



# gradually greater recessive effects (polygenic model)

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


run_with_params 0.001 10.0 0.10 0.25 601
run_with_params 0.001 0.10 0.10 0.25 601
run_with_params 0.001 0.10 1.00 0.25 602
run_with_params 0.001 0.10 10.0 0.25 603
run_with_params 0.001 0.10 99.0 0.25 604

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







