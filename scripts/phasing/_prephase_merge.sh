#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh


readonly input_prefix=${1?Error: Missing arg1 (input_path)}
readonly max_inteval_idx=${2?Error: Missing arg2 (input_path)}
readonly out_merged_file=${3?Error: Missing arg3 (input_path)}



# create list of chunks


# merge chunks


# tabix final merged file





