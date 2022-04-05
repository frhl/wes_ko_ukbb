#!/usr/bin/env bash
#
#$ -N polygenic
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/polygenic.log
#$ -e logs/polygenic.errors.log
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

readonly K=1e-1
readonly h2=3e-1
readonly pi=0
readonly seed=42

readonly out_dir="data/simulation/absence_of_effect_test"
readonly out_prefix="${out_dir}/ukb_eur_h2_${h2}_pi_${pi}_K_${K}_chr${chr}"
readonly out_phenotypes="${out_prefix}_phenotype.tsv.gz"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

simulate_phenotypes() {

  local SECONDS=0
  if [ ! -f "${out_prefix}.tsv.gz" ]; then
    module purge
    set_up_hail
    set_up_pythonpath_legacy
    python3 "${hail_script}" \
       --in_prefix "${in_prefix}"\
       --in_type "mt" \
       --h2 ${h2} \
       --K ${K} \
       --seed ${seed} \
       --simulations 5 \
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


simulate_phenotypes









