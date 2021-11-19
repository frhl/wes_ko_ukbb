Rscript "${rscript}" \
     --out_prefix ${out_prefix}
#!/usr/bin/env bash
#
#$ -N knockout_saige_stats
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/knockout_saige_stats.log
#$ -e logs/knockout_saige_stats.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 8
#$ -q short.qc
#$ -V

#set -o errexit
#set -o nounset

source utils/bash_utils.sh

readonly out_dir="derived/tables/assoc"
readonly out_prefix="${out_dir}/ukb_wes200k"

readonly rscript="scripts/05_knockout_saige_stats.R"

mkdir -p ${out_dir}
set_up_rpy
Rscript "${rscript}" \
    --out_prefix ${out_prefix}


