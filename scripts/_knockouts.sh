#!/usr/bin/env bash
#
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -P lindgren.prjc
#$ -q short.qa

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/05_knockouts.py"
readonly spark_dir="data/tmp/spark"

readonly input_path=${1?Error: Missing arg1 (input_path)}
readonly input_type=${2?Error: Missing arg2 (input_type)}
readonly af_min=${3?Error: Missing arg3 (af_max)}
readonly af_max=${4?Error: Missing arg4 (af_max)}
readonly maf_min=${5?Error: Missing arg5 (maf_min)}
readonly maf_max=${6?Error: Missing arg6 (maf_max)}
readonly exclude=${7?Error: Missing arg7 (Exclude variants)}
readonly in_sex=${8?Error: Missing arg8 (sex)}
readonly in_category=${9?Error: Missing arg9 (variant category)}
readonly only_vcf=${10?Error: Missing arg10 (Should only VCF be created?)}
readonly aggr_method=${11?Error: Missing arg11 (Aggr method: fast or collect)}
readonly out_prefix=${12?Error: Missing arg12 (path prefix for saige output)}
readonly out_type=${13?Error: Missing arg13 (output type e.g., mt,vcf or plink)}

readonly chr=${SGE_TASK_ID}
readonly input_path_chr=$(echo ${input_path} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

evaluate_knockouts() {
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
      --chrom ${chr} \
      --input_path ${input_path_chr} \
      --input_type ${input_type} \
      --csqs_category ${in_category} \
      --use_loftee \
      --export_all_gts \
      ${only_vcf:+--only_vcf} \
      ${aggr_method:+--aggr_method "${aggr_method}"} \
      ${af_min:+--af_min "$af_min"} \
      ${af_max:+--af_min "$af_max"} \
      ${maf_max:+--maf_max "$maf_max"} \
      ${maf_min:+--maf_min "$maf_min"} \
      ${exclude:+--exclude "$exclude"} \
      ${in_sex:+--sex "$in_sex"} \
      --out_prefix ${out_prefix_chr} \
      --out_type ${out_type} \
      && print_update "Finished evaluation knockouts for chr${chr}" ${SECONDS} \
      || raise_error "Evaluating knockouts for chr${chr} failed"
  rm -rf "${out_prefix_chr}_checkpoint.mt"
}


# Find knockouts in data
if [ ! -f "${out_prefix_chr}.vcf.bgz" ]; then
  evaluate_knockouts
else
  >&2 echo "${out_prefix_chr}.vcf.bgz already exists. Skipping.."
fi 

# index resulting knockout VCF for future SAIGE analysis
if [ ! -f "${out_prefix_chr}.vcf.csi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix_chr}.vcf.bgz" "csi"
else
  >&2 echo "${out_prefix_chr}.vcf.bgz.csi (index) already exists. Skipping.."
fi






