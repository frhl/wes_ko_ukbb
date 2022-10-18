#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly input_prefix=${1?Error: Missing arg1 (input_path)}
readonly max_idx=${2?Error: Missing arg2 (input_path)}
readonly out_file=${3?Error: Missing arg3 (input_path)}

readonly mrg_file="${out_file}.tmp"


rm -f ${mrg_file}
# create list of chunks
#for idx in ${1..${max_idx}}; do
for idx in {1..2}; do
  outfile_wo_prefix="${input_prefix}.${idx}of${max_idx}.phased"
  if [ -f "${outfile_wo_prefix}.vcf.gz" ]; then
    echo -e "${outfile_wo_prefix}" >> ${mrg_file}
  else
    raise_error "${outfile_wo_prefix}.vcf.gz does not exist!"
  fi
done





#if [ ! -f ${out_file} ]; then
#  # merge chunks
#  module load BCFtools/1.12-GCC-10.3.0
#  bcftools merge -l ${mrg_file} -Oz -o ${out_file}
#fi

# tabix final merged file
make_tabix "${out_splitted_vcf}" "tbi"





