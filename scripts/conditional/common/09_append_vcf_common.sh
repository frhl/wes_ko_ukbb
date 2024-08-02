#!/usr/bin/env bash
#
# Append pseudovariants and common markers
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=append_vcf_common
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/append_vcf_common.log
#SBATCH --error=logs/append_vcf_common.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 3
#SBATCH --array=1-22

source utils/vcf_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/common/09_append_vcf_common.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly ko_dir="data/knockouts/alt/pp90/combined"
readonly common_dir="data/conditional/common/markers/2024"
readonly out_dir="data/conditional/common/combined/2024"

readonly ko_path_wo_ext="${ko_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense"
readonly common_path_wo_ext="${common_dir}/common_conditional"

readonly ko_path_mt="${ko_path_wo_ext}.mt"
readonly ko_path_vcf="${ko_path_wo_ext}.vcf.bgz"
readonly common_path_mt="${common_path_wo_ext}.mt"
readonly common_path_vcf="${common_path_wo_ext}.vcf.bgz"

readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_w_common"
readonly tmp_markers="${out_prefix}_markers.tmp"
readonly out_markers="${out_prefix}_markers.txt"
readonly out_markers_sorted="${out_prefix}_markers_sorted.txt"

readonly markers_common="${common_dir}/common_conditional.markers"

readonly ko_type="mt"
readonly common_type="mt"
readonly out_type="vcf"
readonly out="${out_prefix}.vcf.bgz"

mkdir -p ${out_dir}

# create file that contains all the markers
if [ ! -f "${out_markers}" ]; then
  cat ${markers_common} | grep "chr${chr}:" > ${tmp_markers}
  if [ "$(cat ${tmp_markers} | wc -l)" -gt "0" ]; then
    mv ${tmp_markers} ${out_markers}
  else
    rm ${tmp_markers}
  fi
fi

# count markers in each file
readonly n_common_markers=$(cat "${out_markers}" | grep chr${chr} | wc -l )
echo "Note: ${n_common_markers} common markers found."

# symlink file if no common markers exist.
if [ "${n_common_markers}" -eq "0" ]; then
  >&2 echo "Note: No common markers. Creating symlink to vcf."
  ln -s "$(pwd)/${ko_path_vcf}" "${out_prefix}.vcf.bgz"
  ln -s "$(pwd)/${ko_path_vcf}.csi" "${out_prefix}.vcf.bgz.csi"
else
  if [ ! -f "${out_prefix}.vcf.bgz" ]; then
    SECONDS=0
    set_up_hail
    set_up_pythonpath_legacy
    python3 "${hail_script}" \
       --chrom ${chr} \
       --ko_path ${ko_path_mt} \
       --ko_type ${ko_type} \
       --common_path ${common_path_mt} \
       --common_type ${common_type} \
       --out_prefix ${out_prefix} \
       --out_type ${out_type} \
       && print_update "Finished merging knockouts with icommon markers ${out_prefix}" ${SECONDS} \
       || raise_error "Merging knockouts with common markers for ${out_prefix} failed!"
  fi

  if [ ! -f "${out_prefix}.vcf.csi" ] & [ "${out_type}" == "vcf" ]; then
    module purge
    module load BCFtools/1.12-GCC-10.3.0
    make_tabix "${out_prefix}.vcf.bgz" "csi"
  fi
fi


