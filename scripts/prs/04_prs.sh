#!/usr/bin/env bash
#
#$ -N prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prs.log
#$ -e logs/prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hge
#$ -t 1
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly rscript="scripts/prs/04_prs.R"
readonly chr=$( get_chr ${SGE_TASK_ID} )

readonly bed_ld_dir="data/prs/hapmap/ld"
readonly bed_ld="${plink_ld_dir}/short_merged_ukb_hapmap_rand_10k_eur.bed"

readonly bed_pred_dir="data/prs/hapmap"
readonly bed_pred="${bed_pred_dir}/ukb_hapmap_500k_eur_chr${chr}.bed"



readonly out_dir="data/prs/test"
readonly out_refix="${out_dir}/long_ukb_ld_hapmap_rand_10k_eur"

mkdir -p ${out_dir}

set_up_rpy
set -x
Rscript "${rscript}" \
   --chrom "AUTOSOMES" \
   --dataset "imp" \
set +x




