#!/usr/bin/env bash
#
#$ -N _merge_grm
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_merge_grm.log
#$ -e logs/_merge_grm.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 6
#$ -q short.qc@@short.hge
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/_merge_grm.py"
readonly spark_dir="data/tmp/spark_dir"

readonly in_prefix=${1?Error: Missing arg2 (in_vcf)}
readonly in_type=${2?Error: Missing arg2 (in_vcf)}
readonly out_prefix=${3?Error: Missing arg2 (in_vcf)}
readonly out_type=${4?Error: Missing arg2 (in_vcf)}

mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy
set -x
python3 "${hail_script}" \
   --in_prefix "${in_prefix}" \
   --in_type "${in_type}" \
   --out_prefix "${out_prefix}" \
   --out_type "${out_type}"
set +x



