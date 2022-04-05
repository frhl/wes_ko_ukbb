#!/usr/bin/env bash
#
#$ -N polygenic_effect
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/polygenic_effect.log
#$ -e logs/polygenic_effect.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hge
#$ -t 21
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/simulation/01_absence_of_effect.py"
readonly spark_dir="data/tmp/spark_dir"

readonly chr=$( get_chr ${SGE_TASK_ID} )

readonly in_dir="data/simulation/mt"
readonly in_prefix="${in_dir}/ukb_eur_25k_samples_chr${chr}.mt"

readonly out_dir="data/simulation/effect_absence"
readonly out_prefix="${out_dir}/ukb_eur_25k_h2_2e-1_pi_None_chr${chr}"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.bed" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --in_prefix "${in_prefix}"\
     --in_type "mt" \
     --chrom "${chr}" \
     --h2 0.2 \
     --seed 42 \
     --simulations 20 \
     --csqs_category "pLoF,damaging_missense" \
     --export_single_markers \
     --out_prefix "${out_prefix}" \
     --out_type "vcf" 
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi




