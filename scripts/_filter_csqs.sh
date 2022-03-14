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

readonly hail_script="scripts/04_filter_csqs.py"
readonly spark_dir="data/tmp/spark"

readonly input_path=${1?Error: Missing arg1 (input_path)}
readonly input_type=${2?Error: Missing arg2 (input_type)}
readonly maf_min=${3?Error: Missing arg3 (maf_min)}
readonly maf_max=${4?Error: Missing arg4 (maf_max)}
readonly in_sex=${5?Error: Missing arg5 (in_sex)}
readonly in_category=${6?Error: Missing arg6 (in_category)}
readonly out_prefix=${7?Error: Missing arg7 (path prefix for saige output)}
readonly out_type=${8?Error: Missing arg8 (output type e.g., mt,vcf or plink)}
readonly exclude=${9?Error: Missing arg13 (Exclude variants)}

readonly chr=${SGE_TASK_ID}
readonly input_path_chr=$(echo ${input_path} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_csq=$(echo "${out_prefix_chr}_${in_category}" | tr "," "_")

echo $out_prefix_csq
echo $out_prefix_chr

if [ ! -f "${out_prefix_csq}.vcf.bgz" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
      --input_path ${input_path_chr} \
      --input_type ${input_type} \
      --csqs_category ${in_category} \
      ${maf_max:+--maf_max "$maf_max"} \
      ${maf_min:+--maf_min "$maf_min"} \
      ${in_sex:+--sex "$in_sex"} \
      ${aggr_method:+--aggr_method "$aggr_method"} \
      ${exclude:+--exclude "$exclude"} \
      --use_loftee \
      --out_prefix ${out_prefix_chr} \
      --out_type ${out_type} \
      && print_update "Finished filtering csqs for chr${chr}" ${SECONDS} \
      || raise_error "Filtering csqs for chr${chr} failed"
else
  >&2 echo "${out_prefix_csq}.vcf.bgz already exists. Skipping.."
fi 




