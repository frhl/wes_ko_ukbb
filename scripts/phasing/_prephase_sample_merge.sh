#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly mergelist=${1?Error: Missing arg3 (input_path)}
readonly expected_files=${2?Error: Missing arg3 (input_path)}
readonly out_file=${3?Error: Missing arg3 (input_path)}

readonly tmplist="${mergelist}.tmp"

cat ${mergelist} | sort | uniq > ${tmplist}
module load BCFtools/1.12-GCC-10.3.0
bcftools merge -l ${tmplist} -Oz -o ${out_file}
rm ${tmplist}




