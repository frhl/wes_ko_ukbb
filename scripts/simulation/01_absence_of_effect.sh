#!/usr/bin/env bash
#
#$ -N absence_of_effect
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/absence_of_effect.log
#$ -e logs/absence_of_effect.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hge
#$ -t 21
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/simulation/simulate_phenotype.py"
readonly spark_dir="data/tmp/spark_dir"

readonly chr=$( get_chr ${SGE_TASK_ID} )

readonly in_dir="data/mt/annotated"
#readonly in_prefix="${in_dir}/ukb_eur_25k_samples_chr${chr}.mt"
readonly in_prefix="${in_dir}/ukb_eur_wes_200k_annot_chr${chr}.mt"

readonly out_dir="data/simulation/effect_absence"
#readonly out_prefix="${out_dir}/ukb_eur_25k_h2_0_pi_None_chr${chr}"
readonly out_prefix="${out_dir}/ukb_eur_h2_0_pi_None_chr${chr}"
readonly out_phenotypes="${out_dir}/ukb_eur_h2_0_pi_None_chr${chr}_phenotype.tsv.gz"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --in_prefix "${in_prefix}"\
     --in_type "mt" \
     --chrom "${chr}" \
     --h2 0 \
     --seed 42 \
     --simulations 20 \
     --csqs_category "pLoF,damaging_missense" \
     --out_prefix "${out_prefix}" \
     --out_type "vcf" \
     && print_update "Finished simulating phenotypes for chr${chr}" ${SECONDS} \
     || raise_error "Simulating phenotypes for for chr${chr} failed"
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi

if [ ! -f "${out_phenotypes}" ]; then
  module purge
  set_up_rpy
  Rscript "${rscript}" \
    --input_path "${out_prefix}.tsv.gz" \
    --real_phenotype_path "${out_prefix}.tsv.gz" \
    --output_path "${out_phenotypes}" \
    && print_update "Finished merging with true phenotypes for chr${chr}" ${SECONDS} \
    || raise_error "Merging with true phenotypes for chr${chr} failed"
fi














