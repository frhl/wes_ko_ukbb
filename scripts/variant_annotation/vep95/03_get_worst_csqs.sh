#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=vep95_get_worst_csqs
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_worst_csqs.log
#SBATCH --error=logs/get_worst_csqs.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1
#SBATCH --dependency="afterok:38762019"

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/vep/vep95/process_csqs"
readonly in="${in_dir}/UKB.chr${chr}.exome_array.variants_only.vep95.csqs.ht"

readonly out_dir="data/vep/vep95/worst_csqs"
readonly out_prefix="${out_dir}/UKB.chr${chr}.exome_array.variants_only.vep95.csqs.worst_csq_by_gene_canonical.original"
readonly hail_script="scripts/variant_annotation/vep95/03_get_worst_csqs.py"

mkdir -p ${out_dir}
mkdir -p ${spark_dir}

set_up_hail 0.2.97
set_up_pythonpath_legacy
python3 ${hail_script} \
     --vep_path "${in}" \
     --out_prefix "${out_prefix}"


