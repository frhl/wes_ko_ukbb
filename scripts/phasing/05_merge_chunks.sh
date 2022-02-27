#!/usr/bin/env bash
#
#$ -N merge_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/merge_chunks.log
#$ -e logs/merge_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qe
#$ -t 20
#$ -V

set -o errexit
set -o nounset
module purge

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly hail_script="scripts/phasing/05_merge_chunks.py"
readonly spark_dir="data/tmp/spark"

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly main_dir="data/phased/wes_union_calls"
readonly in_dir="${main_dir}/ukb_wes_union_calls_200k_chr${chr}-24xshort.qa"
readonly in_prefix="prs52000_pro13000_mprs100000"

readonly out_dir="data/phased/wes_union_calls/test"
readonly out_prefix="merged_test_${chr}"

mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy
SECONDS=0
python3 ${hail_script} \
    --in_dir "${in_dir}" \
    --in_ext ".vcf.gz" \
    --in_prefix "${in_prefix}" \
    --out_prefix "${out_prefix}" \
    --out_type "vcf" \
    && print_update "Finished merging phased data for chr${chr}" ${SECONDS} \
    || raise_error "Merging phased data for chr${chr} failed" 





