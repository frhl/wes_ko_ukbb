#!/usr/bin/env bash
#
# Append pseudo variants with actual variants for downstream conditional analysis.
#
#$ -N combine_ko_rare_common
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/combine_ko_rare_common.log
#$ -e logs/combine_ko_rare_common.errors.log
#$ -P lindgren.prjc
#$ -q short.qc
#$ -pe shmem 3
#$ -t 14-22
#$ -V


set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/combined/01_combine_ko_rare_common.py"

readonly chr="${SGE_TASK_ID}"
readonly variants_dir="data/mt/annotated"

# note: assuming ko and rare variants have already been merged
readonly ko_rare_dir="data/conditional/rare/combined"
readonly common_dir="data/conditional/common/marker_mt"
readonly out_dir="data/conditional/combined"

readonly ko_rare_path_wo_ext="${ko_rare_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense"
readonly common_path_wo_ext="${common_dir}/conditional_markers_chrundefined"

readonly ko_rare_path_mt="${ko_rare_path_wo_ext}.mt"
readonly ko_rare_path_vcf="${ko_rare_path_wo_ext}.vcf.bgz"
readonly common_path_mt="${common_path_wo_ext}.mt"
readonly common_path_vcf="${common_path_wo_ext}.vcf.bgz"

readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense"
readonly out_markers="${out_prefix}_markers.txt"

readonly markers_rare="${ko_rare_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense_markers.txt.gz"
readonly markers_common="${common_dir}/conditional_markers_chrundefined_markers.txt.gz"

readonly ko_rare_type="mt"
readonly common_type="mt"
readonly out_type="vcf"

mkdir -p ${out_dir}

# create new marker file (common markers + rare markers)
zcat $markers_rare | head -n1 > ${out_markers}
zcat $markers_common $markers_rare | grep -v locus | grep chr${chr} >> ${out_markers}
gzip ${out_markers}

# count markers in each file
readonly n_common_markers=$(zcat ${markers_common} | grep chr${chr} | wc -l )
readonly n_rare_markers=$(zcat ${markers_rare} | grep chr${chr} | wc -l )

# if no common markers exists, just copy the previous vcf
# note: we assume that rare markers always exists.
if [ "${n_common_markers}" -eq "0" ]; then
  >&2 echo "Note: No common markers. Creating symlink to rare variants vcf."
  ln -s "$(pwd)/${ko_rare_path_vcf}" "${out_prefix}.vcf.bgz"
  ln -s "$(pwd)/${ko_rare_path_vcf}.csi" "${out_prefix}.vcf.bgz.csi"
elif [ "${n_common_markers}" -eq "0" ]; then
  >&2 echo "Note: No rare markers. Creating symlink to common variants vcf."
  ln -s "$(pwd)/${common_path_vcf}" "${out_prefix}.vcf.bgz"
  ln -s "$(pwd)/${common_path_vcf}.csi" "${out_prefix}.vcf.bgz.csi"
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

