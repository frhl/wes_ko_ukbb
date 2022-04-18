#!/usr/bin/env bash
#
#$ -N _rm_bfile
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_rm_bfile.log
#$ -e logs/_rm_bfile.errors.log
#$ -V

readonly pred=${1?Error: Missing arg3 (prediction file)}
readonly prefix=${2?Error: Missing arg8 (prefix)}

readonly chr="${SGE_TASK_ID}"
readonly out_prefix_chr=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")
readonly tmp_bfile="bfile_${out_prefix_chr}"

echo "Removing backing files.."
rm -f  "${tmp_bfile}.bk"
rm -f "${tmp_bfile}.rds"



