#!/usr/bin/env bash
#
#$ -N merge_markers
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/merge_markers.log
#$ -e logs/merge_markers.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc
#$ -t 3

source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly marker_dir=""
readonly mt_dir=""
readonly out_dir="data/conditional/common"

readonly mt="${mt_dir}/"
readonly markers="${marker_dir}/"

readonly merge_script="scripts/condotional/_merge_markers.sh"

readonly out_prefix="${out_dir}/"

set_up_hail
set_up_pythonpath_legacy
mkdir -p ${out_dir}

submit_merge_markers_job() 
{
  set -x
  qsub -N "_ko_${qsub_name}" \
    -t 6 \
    -q "short.qa" \
    -pe shmem 2 \
    "${merge_script}" \
    "${mt}" \
    "${marker_table}" \
    "${out_prefix}"
  set +x
}


annotation='synonymous'
out_prefix="${out_dir}/211111_intervals_${annotation}_${phenotype}"
merge_markers ${annotation} ${out_prefix}



