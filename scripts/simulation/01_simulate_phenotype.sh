#!/usr/bin/env bash
#
#$ -N simulate_phenotype
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/simulate_phenptype.log
#$ -e logs/simulate_phenptype.errors.log
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

  local K=${1}
  local h2_snp=${2}
  local pi_snp=${3}
  local h2_ch=${4}
  local pi_ch=${5}

  local out_prefix="${out_dir}/ukb_eur_h2_${h2_snp}_${h2_ch}_pi_${pi_snp}_${pi_ch}_K_${K}_chr${chr}"
  local out_phenotypes="${out_prefix}_phenos.tsv.gz"
  
  local sim_name="_sim_${SGE_TASK_ID}"
  local mrg_name="_mrg_${SGE_TASK_ID}"

  set -x
  qsub -N "${sim_name}" \
       -t ${tasks} \
      -q ${queue} \
       -pe shmem ${nslots} \
       ${bash_script} \
       ${in_prefix}\
       ${in_type} \
       ${h2_snp} \
       ${pi_ch} \
       ${h2_snp} \
       ${pi_ch} \
       ${K} \
       ${seed} \
       ${out_prefix}

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
readonly tasks=1-10
readonly seed=42

# simulate traits with no heritability
simulate_phenotypes 1e-1 0 NA NA NA
#simulate_phenotypes 1e-1 1e-1 NA NA NA
simulate_phenotypes 1e-1 1e-3 NA NA NA

simulate_phenotypes 1e-1 0 NA 0 NA
#simulate_phenotypes 1e-1 0 NA 1e-1 NA
simulate_phenotypes 1e-1 0 NA 3e-1 NA

simulate_phenotypes 1e-1 2e-1 NA 3e-1 NA

# simulate traits slightly polygenic traits
#simulate_phenotypes 1e-1 1e-1 0

# simulate moderately polygenic traits
#simulate_phenotypes 1e-1 3e-1 0

# simulate 0.1 % causal varaints
#simulate_phenotypes 1e-1 2e-1 1e-3







