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

readonly hail_script="scripts/simulation/01_absence_of_effect.py"
readonly spark_dir="data/tmp/spark_dir"

readonly chr=$( get_chr ${SGE_TASK_ID} )

readonly in_dir="data/simulation/data"
readonly in_prefix="${in_dir}/ukb_eur_100000_samples_chr${chr}.mt"

readonly out_dir="data/simulation/effect_absence"
readonly out_prefix="${out_dir}/ukb_eur_10k_h2_0_chr${chr}"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.bed" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --h2 "0" \
     --simulations "20" \
     --csqs_category "pLoF,damaging_missense" \
     --in_prefix "${in_prefix}"\
     --in_type "mt" \
     --out_prefix "${out_prefix}" \
     --out_type "vcf" 
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi




