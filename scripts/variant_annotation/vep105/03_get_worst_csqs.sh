#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_worst_csqs
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_worst_csqs.log
#SBATCH --error=logs/get_worst_csqs.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22
#SBATCH --dependency="afterok:38316097"

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/vep/vep105/process_csqs"
readonly in="${in_dir}/UKB.chr${chr}.exome_array.variants_only.vep.csqs.ht"

readonly spliceai_dir="/well/lindgren/barney/brava_annotation/data/spliceai"
readonly spliceai_path="${spliceai_dir}/ukb_wes_450k.spliceai.chr${chr}.ht"

readonly revel_cutoff=0.6
readonly cadd_cutoff=20
readonly case_builder="original" # original/brava
readonly group="original"

readonly out_dir="data/vep/vep105/worst_csqs"
readonly out_prefix="${out_dir}/UKB.chr${chr}.exome_array.variants_only.vep.csqs.worst_csq_by_gene_canonical.${group}"
readonly hail_script="scripts/variant_annotation/vep105/03_get_worst_csqs.py"

mkdir -p ${out_dir}
mkdir -p ${spark_dir}

set_up_hail 0.2.97
set_up_pythonpath_legacy
python3 ${hail_script} \
     --vep_path "${in}" \
     --spliceai_path "${spliceai_path}" \
     --case_builder "${case_builder}" \
     --revel_score "${revel_cutoff}" \
     --cadd_score "${cadd_cutoff}" \
     --out_prefix "${out_prefix}"


