#!/usr/bin/env bash
#
#$ -N filter_plink
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/filter_plink.log
#$ -e logs/filter_plink.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hga
#$ -t 7-8
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/prs/00_filter_plink.py"
readonly spark_dir="data/tmp/spark_dir"

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly input_dir="data/prs/hapmap/ld/old/chr7_and_8"
readonly input_prefix="${input_dir}/short_ukb_hapmap_rand_10k_eur_chr${chr}"
readonly out_dir="data/prs/hapmap/ld"
readonly out_prefix="${out_dir}/short_ukb_hapmap_rand_10k_eur_chr${chr}"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.bed" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --input_prefix ${input_prefix} \
     --input_type "plink" \
     --only_valid_contigs \
     --out_prefix "${out_prefix}" \
     --out_type "plink" 
  set +x
else
  print_update "file ${out_prefix} already exists. Skipping!"
fi




