#!/usr/bin/env bash
#
# @description liftover calls and flip to correct reference sequence.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=liftover_calls
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/liftover_calls.log
#SBATCH --error=logs/liftover_calls.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-21
#
#$ -N liftover_calls
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/liftover_calls.log
#$ -e logs/liftover_calls.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qe
#$ -t 22
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/prefilter/02_liftover_calls.py"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly out_dir="data/unphased/calls/liftover"
readonly out_prefix="${out_dir}/ukb_liftover_calls_500k_chr${chr}"
readonly out_type="mt"

readonly fields_to_drop="fam_id,pat_id,mat_id,is_female,is_case,cm_position"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.mt/_SUCCESS" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --drop_fields ${fields_to_drop} \
     --filter_incorrect_reference \
     --dataset "calls" \
     --liftover
fi





