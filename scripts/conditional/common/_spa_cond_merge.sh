#!/usr/bin/env bash
#
#$ -N _spa_cond_merge
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_spa_cond_merge.log
#$ -e logs/_spa_cond_merge.errors.log
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rmerge="scripts/conditional/common/_spa_cond_merge.R"

readonly prefix=${1?Error: Missing arg1 (prefix)}
readonly out=${2?Error: Missing arg2 (file out full name)}
readonly remove_by_chr=${3?Error: Missing arg3 (Remove indiviudal chromosome files)}

# set up directories 
readonly in_dir=$( dirname ${prefix} )
readonly out_wo_chr=$(echo ${out} | sed -e "s/_chrCHR//g")

# add regex to files we are looking for
readonly file=$(echo ${prefix} | sed -e "s/CHR/[0-9]+/g" )
readonly files="${file##*/}$"

# merge the files using R
set_up_rpy
Rscript ${rmerge} \
  --in_dir ${in_dir} \
  --regex ${files} \
  --out ${out_wo_chr}



