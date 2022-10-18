#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_sample
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/export_sample.log
#SBATCH --error=logs/export_sample.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=21
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/phasing/experimental/04_export_sample.py"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} ) 
readonly in_dir="data/unphased/wes_union_calls"
readonly input_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly input_type="vcf"

readonly eid="1281289,5101274"

readonly out_dir="data/reads/samples"
readonly out_prefix="${out_dir}/eid${eid/","/"_"}_eur_wes_union_calls_chr${chr}"
readonly out_type="vcf"


mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy  
python3 "${hail_script}" \
   --in_file "${input_prefix}" \
   --in_type "${input_type}" \
   --out_prefix ${out_prefix} \
   --out_type "${out_type}" \
   --extract_samples ${eid}



module purge
module load BCFtools/1.12-GCC-10.3.0
make_tabix "${out_prefix}.vcf.bgz" "csi"
make_tabix "${out_prefix}.vcf.bgz" "tbi"

