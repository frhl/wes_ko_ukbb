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
readonly in_dir="data/mt/annotated"
readonly in_prefix="${in_dir}/ukb_eur_wes_200k_annot_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/simulation/phenotypes"

mkdir -p ${out_dir}


simulate_phenotypes() {

  local K=0.1
  local h2_nc=${1}
  local h2_co=${2}
  local h2_ko=${3}
  local pi_nc=${4}
  local pi_co=${5}
  local pi_ko=${6}
  local alpha=${7}
  local beta=${8}
  local theta=${9}

  local h2s="${h2_nc}_${h2_co}_${h2_ko}"
  local pis="${pi_nc}_${pi_co}_${pi_ko}"
  local effects="a${alpha}_b${beta}_t${theta}"

  local out_prefix="${out_dir}/ukb_eur_h2_${h2s}_pi_${pis}_K${K}_${effects}_chr${chr}"
  local out_phenotypes="${out_prefix}_phenos.tsv.gz"
  
  local sim_name="sim_c${SGE_TASK_ID}_h2_${h2s}_pi_${pis}_${effects}"
  local mrg_name="mrg_c${SGE_TASK_ID}_h2_${h2s}_pi_${pis}_${effects}"

  set -x
  if [ 1 -eq 1 ]; then
  qsub -N "${sim_name}" \
       -t ${tasks} \
      -q ${queue} \
       -pe shmem ${nslots} \
       ${bash_script} \
       ${in_prefix}\
       ${in_type} \
       ${h2_nc} \
       ${h2_co} \
       ${h2_ko} \
       ${pi_nc} \
       ${pi_co} \
       ${pi_ko} \
       ${alpha} \
       ${beta} \
       ${theta} \
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
readonly nslots="3"
readonly tasks=1-25
readonly seed=42

# simulate absence of CH effects
simulate_phenotypes 0.00 0.00 0.00 0.00 0.00 0.00 NA NA NA
simulate_phenotypes 0.10 0.00 0.00 0.00 0.00 0.00 NA NA NA
simulate_phenotypes 0.00 0.10 0.00 0.00 0.00 0.00 NA NA NA
simulate_phenotypes 0.10 0.10 0.00 0.00 0.00 0.00 NA NA NA

# simualte CH effect
simulate_phenotypes 0.00 0.00 0.02 0.00 0.00 0.10 NA NA NA
simulate_phenotypes 0.00 0.10 0.02 0.00 0.00 0.10 NA NA NA
simulate_phenotypes 0.00 0.10 0.02 0.00 0.00 0.10 NA NA NA

# simulate effects with betas
simulate_phenotypes 0.00 0.00 0.00 0.00 0.00 0.10 NA NA 0.01
simulate_phenotypes 0.00 0.00 0.00 0.00 0.00 0.10 NA NA 0.10


#simulate_phenotypes 0.10 0.10 0.10 0.10 0.10 0.10 NA NA NA
#simulate_phenotypes 0.00 0.00 0.00 0.10 0.10 0.10 0.01 0.01 0.01
#simulate_phenotypes 0.00 0.00 0.00 0.10 0.10 0.10 0.0 0.0 0.01
#simulate_phenotypes 0.00 0.00 0.00 0.10 0.10 0.10 0.01 0.01 0.00





