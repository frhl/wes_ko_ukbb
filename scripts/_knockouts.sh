#!/usr/bin/env bash
#
# note: this scripts is called from 04_knockouts.sh
#
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_knockouts.log
#$ -e logs/_knockouts.errors.log
#$ -P lindgren.prjc
#$ -q short.qc@@short.hge



set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/04_knockouts.py"
readonly spark_dir="data/tmp/spark"

readonly in_phased=${1?Error: Missing arg1 (in_phased)}
readonly in_phased_type=${2?Error: Missing arg2 (in_phased_type)}
readonly in_unphased=${3?Error: Missing arg3 (in_unphased)}
readonly in_unphased_type=${4?Error: Missing arg4 (in_unphased_type)}
readonly af_max=${5?Error: Missing arg5 (af_max)}
readonly in_category=${6?Error: Missing arg6 (in_category)}
readonly out_prefix=${7?Error: Missing arg7 (path prefix for saige output)}

readonly chr=${SGE_TASK_ID}
readonly phased=$(echo ${in_phased} | sed -e "s/CHR/${chr}/g")
readonly unphased=$(echo ${in_unphased} | sed -e "s/CHR/${chr}/g")
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

set_up_hail
set_up_pythonpath_legacy
set -x
python3 "${hail_script}" \
    --chrom ${chr} \
    --input_phased_path ${phased} \
    --input_unphased_path ${unphased} \
    --input_phased_type ${in_phased_type} \
    --input_unphased_type ${in_unphased_type} \
    --csqs_category ${in_category} \
    --out_prefix ${out} \
    --af_max ${af_max} \
    --export_ko_probability \
    --export_saige_vcf \
    --export_ko_rsid \
    --use_loftee
set +x
print_update "Finished running HAIL for chr${chr}" "${SECONDS}"






