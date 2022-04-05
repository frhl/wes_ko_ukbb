#!/usr/bin/env bash
#
#$ -N simulate_phenotype
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/simulate_phenptype.log
#$ -e logs/simulate_phenptype.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hge
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/bash_utils.sh

readonly hail_script="scripts/simulation/simulate_phenotype.py"
readonly rscript="scripts/simulation/simulate_phenotype.R"
readonly spark_dir="data/tmp/spark_dir"

readonly chr=$( get_chr ${SGE_TASK_ID} )

readonly pheno_dir="data/phenotypes"
readonly pheno_file="${pheno_dir}/filtered_covar_phenotypes_cts.tsv.gz"

readonly in_dir="data/mt/annotated"
readonly in_prefix="${in_dir}/ukb_eur_wes_200k_annot_chr${chr}.mt"

simulate_phenotypes() {

  local K=${1}
  local h2=${2}
  local pi=${3}

  local out_dir="data/simulation/phenotypes"
  local out_prefix="${out_dir}/ukb_eur_h2_${h2}_pi_${pi}_K_${K}_chr${chr}"
  local out_phenotypes="${out_prefix}_phenotype.tsv.gz"

  mkdir -p ${out_dir}
  mkdir -p ${spark_dir}

  local SECONDS=0
  if [ ! -f "${out_prefix}.tsv.gz" ]; then
    echo "Simulating ${simulations} traits (${out_phenotypes})"
    module purge
    set_up_hail
    set_up_pythonpath_legacy
    python3 "${hail_script}" \
       --in_prefix "${in_prefix}"\
       --in_type "mt" \
       --h2 ${h2} \
       --K ${K} \
       --seed ${seed} \
       --simulations ${simulations} \
       --out_prefix "${out_prefix}" \
       && print_update "Finished simulating phenotypes for chr${chr}" ${SECONDS} \
       || raise_error "Simulating phenotypes for for chr${chr} failed"
  else
    print_update "file ${out} already exists. Skipping!"
  fi

  if [ ! -f "${out_phenotypes}" ]; then
    module purge
    set_up_rpy
    Rscript "${rscript}" \
      --input_path "${out_prefix}.tsv.gz" \
      --real_phenotype_path "${pheno_file}" \
      --output_path "${out_phenotypes}" \
      && print_update "Finished merging with true phenotypes for chr${chr}" ${SECONDS} \
      || raise_error "Merging with true phenotypes for chr${chr} failed"
  fi

}

###############
# main script #
###############

readonly seed=42
readonly simulations=30

# simulate traits with no heritability
simulate_phenotypes 1e-1 0 0

# simulate traits slightly polygenic traits
simulate_phenotypes 1e-1 1e-1 0

# simulate moderately polygenic traits
simulate_phenotypes 1e-1 3e-1 0


