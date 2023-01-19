#!/usr/bin/env bash
#
# @description Filter to samples with whole exome sequencing data available. These are the main
# samples that also have knockout status.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=filter_hm3
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/filter_hm3.log
#SBATCH --error=logs/filter_hm3.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 8
#SBATCH --array=1-4
#
#$ -N filter_hm3
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/filter_hm3.log
#$ -e logs/filter_hm3.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -t 1-22
#$ -q long.qc
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/prs/00_filter_hm3.py"
readonly spark_dir="data/tmp/spark_dir"

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly out_dir="data/unphased/imputed/hm3"
readonly out_prefix="${out_dir}/ukb_hm3_500k_chr${chr}"
readonly out_type="mt"

readonly hap_dir="/well/lindgren/flassen/ressources/hapmap/ldpred2"
readonly hap_file="${hap_dir}/map_hm3_plus_liftover.ht"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.mt/_SUCCESS" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --hapmap ${hap_file} \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --chrom "${chr}"
else
  >&2 echo "${out_prefix}.mt already exists. Skipping"
fi



