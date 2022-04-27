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
  local h2_co=${1}
  local pi_co=${2}
  local h2_ko=${3}
  local pi_ko=${4}

  local out_prefix="${out_dir}/ukb_eur_h2_${h2_co}_${h2_ko}_pi_${pi_co}_${pi_ko}_K_${K}_chr${chr}"
  local out_phenotypes="${out_prefix}_phenos.tsv.gz"
  
  local sim_name="sim_${SGE_TASK_ID}"
  local mrg_name="mrg_${SGE_TASK_ID}"


  set -x
  if [ 1 -eq 1 ]; then
  qsub -N "${sim_name}" \
       -t ${tasks} \
      -q ${queue} \
       -pe shmem ${nslots} \
       ${bash_script} \
       ${in_prefix}\
       ${in_type} \
       ${h2_co} \
       ${pi_co} \
       ${h2_ko} \
       ${pi_ko} \
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
readonly tasks=1-50
readonly seed=42

# simulate traits with no heritability
simulate_phenotypes 0.1 0.0 0.09 0.1
#simulate_phenotypes 0.05 0.0 0.05 0.12
#simulate_phenotypes 0.01 0.0 0.1 0.1

#simulate_phenotypes 0.0 0.0 0.0 0.0
#simulate_phenotypes 0.2 0.001 0.0 0.0
#simulate_phenotypes 0.2 0.002 0.0 0.0
#simulate_phenotypes 0.2 0.005 0.0 0.0
#simulate_phenotypes 0.2 0.01 0.0 0.0


# simulate traits slightly polygenic traits
#simulate_phenotypes 1e-1 1e-1 0

# simulate moderately polygenic traits
#simulate_phenotypes 1e-1 3e-1 0

# simulate 0.1 % causal varaints
#simulate_phenotypes 1e-1 2e-1 1e-3







