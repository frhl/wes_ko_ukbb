#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/06_encode_vcf.py"
readonly spark_dir="data/tmp/spark"

readonly input_path=${1?Error: Missing arg1 (input_path)}
readonly input_type=${2?Error: Missing arg2 (input_type)}
readonly in_category=${3?Error: Missing arg3 (variant category)}
readonly only_vcf=${4?Error: Missing arg4 (Should only VCF be created?)}
readonly aggr_method=${5?Error: Missing arg5 (Aggr method: fast or collect)}
readonly out_prefix=${6?Error: Missing arg6 (path prefix for saige output)}
readonly out_type=${7?Error: Missing arg7 (output type e.g., mt,vcf or plink)}

readonly task_id=$( get_array_task_id )
readonly chr=${task_id}

readonly input_path_chr=$(echo ${input_path} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

evaluate_knockouts() {
  echo "Evaluating ${out_prefix_chr}.vcf.bgz.."
  SECONDS=0
  set -x && python3 "${hail_script}" \
      --chrom ${chr} \
      --input_path ${input_path_chr} \
      --input_type ${input_type} \
      --csqs_category ${in_category} \
      --export_all_gts \
      ${only_vcf:+--only_vcf} \
      ${aggr_method:+--aggr_method "${aggr_method}"} \
      --out_prefix ${out_prefix_chr} \
      --out_type ${out_type} \
      --only_singletons \
      --checkpoint
  rm -rf "${out_prefix_chr}_checkpoint.mt"
  rm -rf "${out_prefix_chr}_precheckpoint.mt"
}

set_up_hail
set_up_pythonpath_legacy
if python -c "import hail" &> /dev/null; then
    echo 'Hail can be successfull started.'
else
    >&2 echo 'Error: Hail could not be started.. check conda env?'
    exit 1
fi

evaluate_knockouts

#else
#  >&2 echo "${out_prefix_chr}.vcf.bgz already exists. Skipping!"
#fi


# index resulting knockout VCF for future SAIGE analysis
if [ ! -f "${out_prefix_chr}.vcf.csi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix_chr}.vcf.bgz" "csi"
else
  >&2 echo "${out_prefix_chr}.vcf.bgz.csi (index) already exists. Skipping.."
fi






