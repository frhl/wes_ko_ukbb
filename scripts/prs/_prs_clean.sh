#!/usr/bin/env bash
#
#
#

readonly pred=${1?Error: Missing arg3 (prediction file)}
readonly prefix=${2?Error: Missing arg8 (prefix)}

readonly chr="${SGE_TASK_ID}"
readonly out_prefix_chr=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")
readonly tmp_bfile="${out_prefix_chr}.bfile"

readonly bk="${tmp_bfile}.bk"
readonly rds="${tmp_bfile}.rds"

>&2 echo "Removing ${bk} and ${rds}.."
rm -f ${bk}
rm -f ${rds}





