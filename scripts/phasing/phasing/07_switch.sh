#!/usr/bin/env bash
#
# @description evaluate phasing quality stratified by MAF
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=switch
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/switch.log
#SBATCH --error=logs/switch.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2

#
#$ -N switch
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/switch.log
#$ -e logs/switch.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh

readonly rscript="scripts/phasing/phasing/07_switch.R"

# switch (pp90)
#readonly switch_dir="data/phased/wes_union_calls/200k/shapeit5/switch_pp90"
#readonly switch_regex="ukb_wes_union_calls_200k_shapeit5_chr[0-9]+_pp90.long.mac.new.txt"

# switch (all)
readonly switch_dir="data/phased/wes_union_calls/200k/shapeit5/parents"
readonly switch_regex="ukb_wes_union_calls_200k_shapeit5_parents_chr[0-9]+.long.mac.txt"

# ukb_wes_union_calls_200k_shapeit5_parents_chr8.long.mac.txt.gz

readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/tables/switch_error_rate"
#readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_shapeit5_parents_hail.pp90"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_shapeit5_parents_hail"

readonly wes_variants="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --switch_dir "${switch_dir}" \
  --switch_file "${switch_regex}" \
  --sites "${wes_variants}" \
  --out_prefix "${out_prefix}" 






