#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_worst_csqs
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/get_worst_csqs.log
#SBATCH --error=logs/get_worst_csqs.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-23

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/vep/process_csqs"
readonly in="${in_dir}/UKB.chr${chr}.exome_array.variants_only.vep.csqs.ht"

#readonly alpha_missense_dir="/well/lindgren/flassen/ressources/alphamissense/combined"
#readonly alpha_missense="${alpha_missense_dir}/AlphaMissense.chr${chr}.hg38.ccds.ht"

readonly spliceai_dir="/well/lindgren/barney/brava_annotation/data/spliceai"
readonly spliceai_path="${spliceai_dir}/ukb_wes_450k.spliceai.chr${chr}.ht"

#readonly am_cutoff=0.564
readonly revel_cutoff=0.773
readonly spliceai_cutoff=0.50 # 0.20/0.50
readonly case_builder="brava" # original/brava
readonly group="brava_s50"

readonly out_dir="data/vep/worst_csqs"
readonly out_prefix="${out_dir}/UKB.chr${chr}.exome_array.variants_only.vep.csqs.worst_csq_by_gene_canonical.${group}"
readonly hail_script="scripts/variant_annotation/03_get_worst_csqs.py"

mkdir -p ${out_dir}
mkdir -p ${spark_dir}

set_up_hail 0.2.97
set_up_pythonpath_legacy
python3 ${hail_script} \
     --vep_path "${in}" \
     --case_builder "${case_builder}" \
     --spliceai_path "${spliceai_path}" \
     --spliceai_score "${spliceai_cutoff}" \
     --revel_score "${revel_cutoff}" \
     --out_prefix "${out_prefix}"


