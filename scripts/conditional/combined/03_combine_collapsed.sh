#!/usr/bin/env bash
#
# Append pseudo variants with actual variants for downstream conditional analysis.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=combine_collapsed
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/combine_collapsed.log
#SBATCH --error=logs/combine_collapsed.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/combined/03_combine_collapsed.py"

readonly array_idx=$( get_array_task_id )
readonly chr="$( get_chr ${array_idx} )"
readonly variants_dir="data/mt/annotated"

readonly ko_dir="data/knockouts/alt/pp90/combined"
readonly collapsed_dir="data/mt/dosages_urv/pp90"
readonly rare_dir="data/mt/dosages_urv/pp90"
readonly common_dir="data/conditional/common/markers"

readonly collapsed_path_wo_ext="${collapsed_dir}/ukb_eur_wes_200k_chr${chr}_max_ds"
readonly rare_path_wo_ext="${rare_dir}/ukb_eur_wes_200k_chr${chr}_mac_gt10"
readonly common_path_wo_ext="${common_dir}/common_conditional"

readonly ko_path="${ko_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense.mt"
readonly ko_type="mt"

readonly collapsed_path="${collapsed_path_wo_ext}.mt"
readonly collapsed_type="mt"

readonly rare_path="${rare_path_wo_ext}.mt"
readonly rare_type="mt"

readonly common_path_vcf="${common_path_wo_ext}.vcf.bgz"
readonly common_path="${common_path_wo_ext}.mt"
readonly common_type="mt"

readonly markers_common="${common_path_wo_ext}.markers"

readonly out_dir="data/conditional/combined/combine_collapsed_urv"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense"
readonly out_type="vcf"
readonly out="${out_prefix}.vcf.bgz"

mkdir -p ${out_dir}

# A function to extract all the (unique) chromsomes
extract_chr_from_vcf() {
  >&2 echo "Extracting chromsomes from VCF ${1}.."
  module load BCFtools/1.12-GCC-10.3.0
  local string=$(bcftools query -f '%CHROM\n' ${1} | sort | uniq)
  echo ${string}
}

# check if the a common chromosome is here, only returning exact match
# note to self: need to ensure that this is converted into a bash array
readonly chr_common=$(extract_chr_from_vcf ${common_path_vcf})
readonly chr_common_present=$(echo ${chr_common} | grep -ow "chr${chr}" | wc -l)

if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  set +u
  set_up_hail
  set -u
  set_up_pythonpath_legacy
  if [ "${chr_common_present}" -eq "0" ]; then
    echo "No common markers for chr${chr} in ${common_path_vcf}. Using rare and KOs."
    python3 "${hail_script}" \
       --chrom ${chr} \
       --ko_path ${ko_path} \
       --ko_type ${ko_type} \
       --collapsed_path ${collapsed_path} \
       --collapsed_type ${collapsed_type} \
       --rare_path ${rare_path} \
       --rare_type ${rare_type} \
       --out_prefix ${out_prefix} \
       --out_type ${out_type}
  else
    echo "Found common markers for chr${chr} in ${common_path_vcf}!"
    python3 "${hail_script}" \
       --chrom ${chr} \
       --ko_path ${ko_path} \
       --ko_type ${ko_type} \
       --collapsed_path ${collapsed_path} \
       --collapsed_type ${collapsed_type} \
       --rare_path ${rare_path} \
       --rare_type ${rare_type} \
       --common_path ${common_path} \
       --common_type ${common_type} \
       --out_prefix ${out_prefix} \
       --out_type ${out_type}
  fi
fi

# index the resulting file.
if [ ! -f "${out_prefix}.vcf.csi" ] & [ "${out_type}" == "vcf" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix}.vcf.bgz" "csi"
fi



