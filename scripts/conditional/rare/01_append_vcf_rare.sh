#!/usr/bin/env bash
#
# @description Append pseudo variants with actual variants for downstream conditional analysis.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=append_vcf_rare
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/append_vcf_rare.log
#SBATCH --error=logs/append_vcf_rare.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22


set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/rare/01_append_vcf_rare.py"

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly variants_dir="data/mt/prefilter/pp90"
readonly ko_dir="data/knockouts/alt/pp90/combined"
readonly out_dir="data/conditional/rare/combined/mt"

readonly variants_path="${variants_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.mt"
readonly input_path="${ko_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense.mt"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense"
readonly variants_type="mt"
readonly input_type="mt"
readonly out_type="vcf"

readonly category="pLoF,damaging_missense"
readonly exclude="data/genes/220310_common_plofs_to_exclude.txt"

mkdir -p ${out_dir}

# Note: assuming variants_path is a MatrixTable

if [  -d "${variants_path}" ]; then
  if [ ! -f "${out_prefix}.vcf.bgz" ]; then
    SECONDS=0
    set_up_hail
    set_up_pythonpath_legacy
    python3 "${hail_script}" \
       --ko_path ${input_path} \
       --ko_type ${input_type} \
       --var_path ${variants_path} \
       --var_type ${variants_type} \
       --out_type ${out_type} \
       --out_prefix ${out_prefix} \
       --csqs_category ${category} \
       --exclude $exclude \
       && print_update "Finished merging knockouts with markers ${out_prefix}" ${SECONDS} \
       || raise_error "Merging knockouts with markers for ${out_prefix} failed!"
  fi
else
  >&2 echo "${variants_path} does not exist!"
fi

if [ ! -f "${out_prefix}.vcf.csi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix}.vcf.bgz" "csi"
fi



