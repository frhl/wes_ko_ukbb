#!/usr/bin/env bash
#
#$ -N merge_bed
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/merge_bed.log
#$ -e logs/merge_bed.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 6
#$ -q short.qc@@short.hge
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/prs/hapmap/ld/unrel_kin_eur_10k"
readonly in_prefix="${in_dir}/short_ukb_hapmap_rand_10k_eur_chr"

readonly out_dir="data/prs/hapmap/ld/unrel_kin_eur_10k"
readonly out_prefix="${out_dir}/short_merged_ukb_hapmap_rand_10k_eur"

readonly hail_script="scripts/prs/02_merge_bed.py"
readonly spark_dir="data/tmp/spark_dir"

mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.bed" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --in_prefix "${in_prefix}" \
     --in_type "plink" \
     --out_prefix "${out_prefix}" \
     --out_type "plink" 
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi




