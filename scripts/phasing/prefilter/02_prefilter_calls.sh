#!/usr/bin/env bash
#
# @description generate files of genotyped calls
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prefilter_calls
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prefilter_calls.log
#SBATCH --error=logs/prefilter_calls.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 3
#SBATCH --array=1-22
#
#$ -N prefilter_calls
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prefilter_calls.log
#$ -e logs/prefilter_calls.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qe
#$ -t 21-22
#$ -V


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/prefilter/prefilter.py"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/unphased/calls/liftover"
readonly in_file="${in_dir}/ukb_liftover_calls_500k_chr${chr}.mt"
readonly in_type="mt"

#readonly out_dir_500k="data/unphased/calls/prefilter_no_maf_cutoff/500k"
#readonly out_prefix_500k="${out_dir_500k}/ukb_prefilter_calls_500k_chr${chr}"
#readonly out_type_500k="mt"

# Note: 200k samples that overlap samples in WES
readonly samples_list="data/unphased/overlap/ukb_calls_wes_samples.txt"
readonly out_dir_200k="data/unphased/calls/prefilter_no_maf_cutoff/200k"
readonly out_prefix_200k="${out_dir_200k}/ukb_prefilter_calls_200k_chr${chr}"
readonly out_type_200k="mt"

readonly missing=0.05
readonly ancestry=""

mkdir -p ${spark_dir}
mkdir -p ${out_dir_200k}
mkdir -p ${out_dir_500k}

set_up_hail
set_up_pythonpath_legacy

if [ ! -f "${out_prefix_200k}.mt/_SUCCESS" ]; then
  python3 ${hail_script} \
     --input_path ${in_file} \
     --input_type ${in_type} \
     --out_prefix ${out_prefix_200k} \
     --out_type ${out_type_200k} \
     --extract_samples ${samples_list} \
     --missing ${missing}
fi





