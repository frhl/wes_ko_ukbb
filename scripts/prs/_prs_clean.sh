#!/usr/bin/env bash

source utils/qsub_utils.sh

readonly pred=${1?Error: Missing arg3 (prediction file)}
readonly prefix=${2?Error: Missing arg8 (prefix)}

readonly cluster=$( get_current_cluster )
readonly chr=$( get_array_task_id )

readonly out_prefix_chr=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")
readonly tmp_bfile="${out_prefix_chr}.bfile"

readonly bk="${tmp_bfile}.bk"
readonly rds="${tmp_bfile}.rds"

>&2 echo "Removing ${bk} and ${rds}.."
rm -f ${bk}
rm -f ${rds}





