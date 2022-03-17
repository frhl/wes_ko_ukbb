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
readonly out_prefix=${3?Error: Missing arg7 (path prefix for saige output)}
readonly out_type=${4?Error: Missing arg8 (output type e.g., mt,vcf or plink)}
readonly aggr_method=${5?Error: Missing arg9 (Aggr method: fast or collect)}
readonly randomize_phase=${6?Error: Missing arg10 (Should phase be randomized?)}
readonly seed=${7?Error: Missing arg11 (Seed for random operations)}
readonly only_vcf=${8?Error: Missing arg12 (Only return VCF)}

readonly chr=${SGE_TASK_ID}
readonly input_path_chr=$(echo ${input_path} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

if [ ! -f "${out_prefix_chr}.vcf.bgz" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
      --chrom ${chr} \
      --input_path ${input_path_chr} \
      --input_type ${input_type} \
      ${seed:+--seed "$seed"} \
      ${randomize_phase:+--randomize_phase} \
      ${only_vcf:+--only_vcf} \
      ${aggr_method:+--aggr_method "$aggr_method"} \
      --out_prefix ${out_prefix_chr} \
      --out_type ${out_type} \
      && print_update "Finished calculating knockouts for chr${chr}" ${SECONDS} \
      || raise_error "Calculating knockouts for chr${chr} failed"
else
  >&2 echo "${out_prefix_chr}.vcf.bgz already exists. Skipping.."
fi 

if [ ! -f "${out_prefix_chr}.vcf.csi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix_chr}.vcf.bgz" "csi"
else
  >&2 echo "${out_prefix_chr}.vcf.bgz.csi (index) already exists. Skipping.."
fi






