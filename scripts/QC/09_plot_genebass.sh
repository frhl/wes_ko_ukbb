#!/usr/bin/env bash
#
#$ -N plot_genebass
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/plot_genebass.log
#$ -e logs/plot_genebass.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hge


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="data/variants"
readonly out_dir="derived/variants"

# hail script
readonly rcode="scripts/QC/09_plot_genebass.R"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_genebass_AF_noNA"

# run 
set_up_RSAIGE
mkdir -p ${out_dir}
Rscript "${rcode}" \
     --in_dir ${in_dir} \
     --in_pattern 'genebass' \
     --out_prefix ${out_prefix}


