#!/usr/bin/env bash
#
#$ -N wes_union_calls_gen
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_union_calls_gen.log
#$ -e logs/wes_union_calls_gen.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hga
#$ -t 1-22
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/02_geno_gen.py"
readonly in_dir="data/unphased/wes/post-qc"

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/unphased/wes_union_calls/test"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly out_type="vcf"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}
if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --input_path "${in_file}" \
     --input_type "${in_type}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --min_mac 2 \
     --missing 0.05 \
     --dataset "calls" \
     --liftover
  set +x
  log_runtime ${SECONDS}
else
  print_update "file ${out} already exists. Skipping!"
fi

module purge
module load BCFtools/1.12-GCC-10.3.0
make_tabix "${out_prefix}.vcf.bgz" "tbi"



