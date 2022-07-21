#!/usr/bin/env bash
#
#$ -N simulate_phenotype
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/simulate_phenotype.log
#$ -e logs/simulate_phenotype.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 21
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
readonly in_dir="data/simulation/data"
#readonly in_prefix="${in_dir}/ukb_eur_100k_chr${chr}.mt"
readonly in_prefix="${in_dir}/ukb_eur_10000_samples_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/simulation/phenotypes"

mkdir -p ${out_dir}


simulate_phenotypes() {

  local K=0.1
  local h2=${1}
  local var_beta=${2}
  local var_theta=${3}
  local pi_beta=${4}
  local pi_theta=${5}

  local vars="${var_beta}_${var_theta}"
  local pis="${pi_beta}_${pi_theta}"

  local out_prefix="${out_dir}/ukb_eur_h2_${h2}_var_${vars}_pi_${pis}_K${K}_chr${chr}"
  local out_phenotypes="${out_prefix}_phenos.tsv.gz"
  
  local sim_name="sim_c${SGE_TASK_ID}_h2_${h2}_var_${vars}_pi_${pis}"
  local mrg_name="mrg_c${SGE_TASK_ID}_h2_${h2}_var_${vars}_pi_${pis}"

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
       ${pi_beta} \
       ${pi_theta} \
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
readonly tasks=1-20
readonly seed=42

# simulate absence of CH effects
simulate_phenotypes 0.00 0.00 0.00 0.01 0.01

# gradually greater recessive effects (polygenic model)
simulate_phenotypes 0.30 0.01 0.01 1.00 1.00
simulate_phenotypes 0.30 0.01 0.10 1.00 1.00
simulate_phenotypes 0.30 0.01 0.50 1.00 1.00
simulate_phenotypes 0.30 0.01 1.00 1.00 1.00
simulate_phenotypes 0.30 0.01 10.0 1.00 1.00
simulate_phenotypes 0.30 0.01 20.0 1.00 1.00
simulate_phenotypes 0.30 0.01 40.0 1.00 1.00
simulate_phenotypes 0.30 0.01 100.0 1.00 1.00
simulate_phenotypes 0.30 0.01 500.0 1.00 1.00
simulate_phenotypes 0.30 0.01 1000.0 1.00 1.00





