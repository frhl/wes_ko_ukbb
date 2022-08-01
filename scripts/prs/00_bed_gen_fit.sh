#!/usr/bin/env bash
#
#$ -N bed_gen_fit
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/bed_gen_fit.log
#$ -e logs/bed_gen_fit.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q long.qc@@long.hge
#$ -t 1-21
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/prs/00_bed_gen.py"
readonly spark_dir="data/tmp/spark_dir"
readonly hap_dir="/well/lindgren/flassen/ressources/hapmap"

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly out_dir="data/prs/hapmap/ukb_500k/fitting"
readonly out_prefix="${out_dir}/ukb_hapmap_500k_eur_chr${chr}"
readonly hap_file="${hap_dir}/weights.l2.ldscore.liftover.ht"

readonly sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.bed" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --dataset "imp" \
     --ancestry "eur" \
     --hapmap ${hap_file} \
     --filter_missing 0.01 \
     --exclude_related \
     --exclude_samples "${sample_list}" \
     --filter_to_unrelated_using_kinship_coef \
     --min_maf 0.01 \
     --liftover \
     --dbsnp \
     --out_prefix "${out_prefix}" \
     --out_type "plink" 
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi




