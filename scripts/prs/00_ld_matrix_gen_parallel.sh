#!/usr/bin/env bash
#
#$ -N ld_matrix_gen_parallel
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/ld_matrix_gen_parallel.log
#$ -e logs/ld_matrix_gen_parallel.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hge
#$ -t 1-22
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/prs/00_bed_gen.py"
readonly spark_dir="data/tmp/spark_dir"
readonly hap_dir="/well/lindgren/flassen/ressources/hapmap"

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly out_dir="data/prs/hapmap/ld"
readonly out_prefix="${out_dir}/short_ukb_hapmap_rand_10k_eur_chr${chr}"
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
     --extract_samples "${sample_list}" \
     --random_samples 10000 \
     --hapmap "${hap_file}" \
     --min_maf 0.005 \
     --liftover \
     --dbsnp \
     --out_prefix "${out_prefix}" \
     --out_type "plink" 
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi




