#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=filter_maf_pp
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/filter_maf_pp.log
#SBATCH --error=logs/filter_maf_pp.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22

# this script used to be called 02_prefilter_variants.*

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/02_filter_maf_pp.py"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

#readonly in_dir="data/mt/annotated/new" #(samples=)
readonly in_dir="data/mt/annotated/old" #(samples=176587)
readonly input_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.mt"
readonly input_type="mt"

#readonly out_dir="data/mt/prefilter/pp90"
readonly out_dir="data/mt/prefilter/no_pp_cutoff/old"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005"
readonly out="${out_prefix}.vcf.bgz"
readonly out_type="vcf"

# remove these common plofs (90% pop)
readonly exclude="data/genes/220310_common_plofs_to_exclude.txt"

readonly maf_min=0.00
readonly maf_max=0.05
#readonly pp_cutoff=0.90
readonly partitions=64

mkdir -p ${out_dir}

module purge && module load BCFtools
make_tabix "${out}" "csi"
exit 0

plink2

#if [ ! -f "${out_prefix}.mt/_SUCCESS" ]; then
if [ ! -f "${out}" ]; then
  #rm -rf "${out_prefix}.mt/"
  set_up_hail 0.2.97
  set_up_pythonpath_legacy
  
  SECONDS=0
  python3 "${hail_script}" \
     --input_path ${input_prefix}\
     --input_type ${input_type} \
     --out_prefix ${out_prefix} \
     --out_type ${out_type} \
     --maf_min ${maf_min} \
     --maf_max ${maf_max} \
     --exclude ${exclude} \
     --partitions ${partitions} \
     --export_csqs \
     && print_update "Finished annotating MatrixTables chr${chr}" ${SECONDS} \
     || raise_error "Annotating MatrixTables for chr${chr} failed"
  # ensure that VCF is indexed
  # make_tabix "${out_prefix}.vcf.bgz" "csi"
  # --pp_cutoff ${pp_cutoff} \
  module purge && module load BCFtools
  make_tabix "${out}" "csi"
else
  >&2 echo "${out_prefix}.mt already exists! Skipping"
fi










