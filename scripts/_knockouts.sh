#!/usr/bin/env bash
#
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_knockouts.log
#$ -e logs/_knockouts.errors.log
#$ -P lindgren.prjc
#$ -q short.qa

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/04_knockouts.py"
readonly spark_dir="data/tmp/spark"

readonly input_path=${1?Error: Missing arg1 (input_path)}
readonly input_type=${2?Error: Missing arg2 (input_type)}
readonly af_min=${3?Error: Missing arg3 (af_max)}
readonly af_max=${4?Error: Missing arg4 (af_max)}
readonly maf_min=${5?Error: Missing arg5 (maf_min)}
readonly maf_max=${6?Error: Missing arg6 (maf_max)}
readonly in_sex=${7?Error: Missing arg7 (in_sex)}
readonly in_category=${8?Error: Missing arg8 (in_category)}
readonly out_prefix=${9?Error: Missing arg9 (path prefix for saige output)}
readonly out_type=${10?Error: Missing arg10 (output type e.g., mt,vcf or plink)}
readonly randomize_phase=${11?Error: Missing arg11 (Should phase be randomized? For creating Null knockout model)}

readonly chr=${SGE_TASK_ID}
readonly input_path_chr=$(echo ${input_path} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

SECONDS=0
set_up_hail
set_up_pythonpath_legacy
set -x
python3 "${hail_script}" \
    --chrom ${chr} \
    --input_path ${input_path_chr} \
    --input_type ${input_type} \
    --csqs_category ${in_category} \
    ${af_min:+--af_min "$af_min"} \
    ${af_max:+--af_min "$af_max"} \
    ${maf_max:+--maf_max "$maf_max"} \
    ${maf_min:+--maf_min "$maf_min"} \
    ${in_sex:+--sex "$in_sex"} \
    ${randomize_phase:+--randomize_phase} \
    --use_loftee \
    --out_prefix ${out_prefix_chr} \
    --out_type ${out_type} \
    && print_update "Finished calculating knockouts for chr${chr}" ${SECONDS} \
    || raise_error "Calculating knockouts for chr${chr} failed"

if [ ! -f "${out_prefix_chr}.vcf.tbi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix_chr}.vcf.bgz" "tbi"
fi







