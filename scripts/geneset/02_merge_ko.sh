#!/usr/bin/env bash
#
#$ -N merge_ko
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/merge_ko.log
#$ -e logs/merge_ko.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 6
#$ -q short.qc@@short.hge
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/knockouts/alt"
readonly in_prefix="${in_dir}/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense.mt"

readonly out_dir="data/geneset/knockouts"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_maf0to5e-2_pLoF_damaging_missense_combined"

readonly hail_script="scripts/geneset/02_merge_ko.py"
readonly spark_dir="data/tmp/spark_dir"

mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.bed" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --in_prefix "${in_prefix}" \
     --in_type "mt" \
     --out_prefix "${out_prefix}" \
     --out_type "vcf" 
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi




