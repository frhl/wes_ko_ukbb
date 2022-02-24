#!/usr/bin/env bash
#
#$ -N calls_gen
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/calls_gen.log
#$ -e logs/calls_gen.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 21
#$ -V

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/01_geno_gen.py"

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly out_dir="data/unphased/calls"
readonly out_prefix="${out_dir}/ukb_prefilter_calls_500k_chr${chr}"
readonly out_type="vcf"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "$( get_hail_ext ${out_prefix} ${out_type})" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --liftover \
     --min_mac 2 \
     --missing 0.05 \
     --ancestry "eur" \
     --dataset "calls"
  set +x
  log_runtime ${SECONDS}
else
  echo "file ${out} already exists. Skipping!"
fi


if [[ ${out_type} == "vcf" ]]; then
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix}.vcf.bgz" "tbi"
fi



