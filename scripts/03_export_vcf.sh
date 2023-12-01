#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_vcf
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/export_vcf.log
#SBATCH --error=logs/export_vcf.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-21

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/03_export_vcf.py"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/mt/prefilter/pp90"
readonly input_path="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.vcf.bgz"
readonly input_type="vcf.bgz"

readonly out_dir="data/mt/prefilter/pp90"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt"
readonly out_type="vcf"

mkdir -p ${out_dir}

module load BCFtools
bcftools view -c 1 -c 1:nonmajor ${input_path} | bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' -Oz -o ${out_prefix}.vcf.gz

# deprecated in favour of BCFtools to annotate variants
#set_up_hail
#set_up_pythonpath_legacy
#python3 "${hail_script}" \
#   --input_path ${input_path}\
#   --input_type ${input_type} \
#   --out_prefix ${out_prefix} \
#   --out_type ${out_type} \


 







