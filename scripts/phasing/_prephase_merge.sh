#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly hail_script="scripts/phasing/_prephase_merge.py"

readonly input_list=${1?Error: Missing arg1 (input_list)}
readonly input_type=${2?Error: Missing arg3 (input_type)}
readonly output_prefix=${3?Error: Missing arg3 (output_file)}
readonly output_type=${4?Error: Missing arg4 (output_type)}

readonly tmp="${output_prefix}.tmp"

#rm -f ${mrg_list}
# create list of chunks
#for idx in {1..2}; do
#for idx in {1..${max_idx}}; do
#  outfile_w_prefix="${input_prefix}.${idx}of${max_idx}.phased.vcf.gz"
#  if [ -f "${outfile_w_prefix}" ]; then
#    echo -e "${outfile_w_prefix}" >> ${mrg_list}
#  else
#    raise_error "${outfile_w_prefix} does not exist!"
#  fi
#done


# combine all files into a single one
#module purge
#set_up_hail
#set_up_pythonpath_legacy
#python3 ${hail_script} \
#  --input_list ${input_list} \
#  --input_type ${input_type} \
#  --out_prefix ${output_prefix} \
#  --out_type ${output_type}

# remove duplicates (from debugging the functions)
cat ${input_list} | sort | uniq  > ${tmp}
# combine VCFs
module load BCFtools/1.12-GCC-10.3.0
bcftools merge -l ${tmp} -o "${output_prefix}.vcf"
bgzip "${output_prefix}.vcf"

rm ${tmp}

if [ "${out_type}" == "vcf" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_file}.vcf.gz" "tbi"
fi




