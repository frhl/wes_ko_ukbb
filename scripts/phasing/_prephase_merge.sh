#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly hail_script="scripts/phasing/_prephase_merge.py"

readonly input_prefix=${1?Error: Missing arg1 (input_path)}
readonly max_idx=${2?Error: Missing arg2 (input_path)}
readonly out_file=${3?Error: Missing arg3 (input_path)}

readonly mrg_list="${out_file}.files"
readonly mrg_type="vcf"
readonly out_type="vcf"

rm -f ${mrg_list}
# create list of chunks
#for idx in {1..2}; do
for idx in {1..${max_idx}}; do
  outfile_w_prefix="${input_prefix}.${idx}of${max_idx}.phased.vcf.gz"
  if [ -f "${outfile_w_prefix}" ]; then
    echo -e "${outfile_w_prefix}" >> ${mrg_list}
  else
    raise_error "${outfile_w_prefix} does not exist!"
  fi
done


# combine all files into a single one
module purge
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --input_list ${mrg_list} \
  --input_type ${mrg_type} \
  --out_prefix ${out_file} \
  --out_type ${out_type}

# index resulting vcf
module purge
module load BCFtools/1.12-GCC-10.3.0
make_tabix "${out_file}.vcf.bgz" "tbi"





