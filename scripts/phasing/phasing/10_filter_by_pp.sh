#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=filter_by_pp
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/filter_by_pp.log
#SBATCH --error=logs/filter_by_pp.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-21

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/phasing/phasing/10_filter_by_pp.py"

readonly input_dir="data/phased/wes_union_calls/200k/shapeit5/clean_ligated"
readonly input_path="${input_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly input_type="vcf"

readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/filter_by_pp"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.clean"
readonly out="${out_prefix}.vcf.bgz"
readonly out_type="vcf"

# remove these common plofs (90% pop)
readonly exclude="data/genes/220310_common_plofs_to_exclude.txt"

readonly maf_min=0.00
readonly maf_max=0.05
readonly pp_cutoff=0.90

if [ ! -f "${out}" ]; then
  SECONDS=0
  set_up_hail 0.2.97
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --input_path ${input_path}\
     --input_type ${input_type} \
     --out_prefix ${out_prefix} \
     --out_type ${out_type} \
     --maf_min ${maf_min} \
     --maf_max ${maf_max} \
     --exclude ${exclude} \
     --pp_cutoff ${pp_cutoff} \
     && print_update "Finished annotating MatrixTables chr${chr}" ${SECONDS} \
     || raise_error "Annotating MatrixTables for chr${chr} failed"
  module purge && module load BCFtools
  make_tabix "${out}" "csi"
else
  >&2 echo "${out_prefix}.mt already exists! Skipping"
fi






