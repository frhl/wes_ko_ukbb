#!/usr/bin/env bash
#
# Append pseudo variants with actual variants for downstream conditional analysis.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=combine_collapsed
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/combine_collapsed.log
#SBATCH --error=logs/combine_ko_rare_cond_common.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 5
#SBATCH --array=22
#
#$ -N combine_collapsed
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/combine_collapsed.log
#$ -e logs/combine_collapsed.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc
#$ -t 7-22
#$ -V

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/combined/01_combine_collpased.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )
readonly variants_dir="data/mt/annotated"

# note: assuming ko and rare variants have already been merged
readonly ko_rare_dir="data/conditional/rare/combined/mt"
readonly common_dir="data/conditional/common/markers"

readonly ko_rare_path_wo_ext="${ko_rare_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense"
readonly common_path_wo_ext="${common_dir}/common_conditional"

readonly ko_rare_path_mt="${ko_rare_path_wo_ext}.mt"
readonly ko_rare_path_vcf="${ko_rare_path_wo_ext}.vcf.bgz"
readonly common_path_mt="${common_path_wo_ext}.mt"
readonly common_path_vcf="${common_path_wo_ext}.vcf.bgz"
readonly markers_common="${common_path_wo_ext}.markers"

# one file of rare markers for each chromosome
readonly markers_rare="${ko_rare_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_markers.txt.gz"

readonly out_dir="data/conditional/combined/combine_collpased"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense"

readonly ko_rare_type="mt"
readonly common_type="mt"
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

if [ "${chr_common_present}" -eq "0" ]; then
  echo "No common markers for ${chr} in ${common_path_vcf}. Creating symlink to rare variants"
  ln -s "$(pwd)/${ko_rare_path_vcf}" "${out_prefix}.vcf.bgz"
  ln -s "$(pwd)/${ko_rare_path_vcf}.csi" "${out_prefix}.vcf.bgz.csi"
else 
  if [ ! -f "${out_prefix}.vcf.bgz" ]; then
    SECONDS=0
    set_up_hail
    set_up_pythonpath_legacy
    python3 "${hail_script}" \
       --chrom ${chr} \
       --ko_rare_path ${ko_rare_path_mt} \
       --ko_rare_type ${ko_rare_type} \
       --common_path ${common_path_mt} \
       --common_type ${common_type} \
       --out_prefix ${out_prefix} \
       --out_type ${out_type} \
       && print_update "Finished merging knockouts with markers ${out_prefix}" ${SECONDS} \
       || raise_error "Merging knockouts with markers for ${out_prefix} failed!"
  fi
  # index the resulting file.
  if [ ! -f "${out_prefix}.vcf.csi" ] & [ "${out_type}" == "vcf" ]; then
    module purge
    module load BCFtools/1.12-GCC-10.3.0
    make_tabix "${out_prefix}.vcf.bgz" "csi"
  fi
fi



