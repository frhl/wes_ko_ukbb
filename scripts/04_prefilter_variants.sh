#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prefilter_variants
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prefilter_variants.log
#SBATCH --error=logs/prefilter_variants.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 2
#SBATCH --array=1-22
#SBATCH --requeue
#
#
#$ -N prefilter_variants
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prefilter_variants.log
#$ -e logs/prefilter_variants.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 1-19
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/04_prefilter_variants.py"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/mt/annotated"
readonly input_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.mt"
readonly input_type="mt"

readonly out_dir="data/mt/prefilter"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp95.maf0_005"
readonly out_type="mt"

# remove these common plofs (90% pop)
readonly exclude="data/genes/220310_common_plofs_to_exclude.txt"

readonly maf_min=0.00
readonly maf_max=0.05
readonly pp_cutoff=0.95

mkdir -p ${out_dir}

SECONDS=0
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --input_path ${input_prefix}\
   --input_type ${input_type} \
   --out_prefix ${out_prefix} \
   --out_type ${out_type} \
   --pp_cutoff ${pp_cutoff} \
   --maf_min ${maf_min} \
   --maf_max ${maf_max} \
   --exclude ${exclude} \
   && print_update "Finished annotating MatrixTables chr${chr}" ${SECONDS} \
   || raise_error "Annotating MatrixTables for chr${chr} failed"





