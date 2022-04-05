#!/usr/bin/env bash
#
#$ -N bed_gen
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/bed_gen.log
#$ -e logs/bed_gen.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qc@@short.hge
#$ -t 20
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/simulation/00_bed_gen.py"
readonly spark_dir="data/tmp/spark_dir"

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly in_dir="data/mt/annotated"
readonly in_prefix="${in_dir}/ukb_eur_wes_200k_annot_chr${chr}.mt"

readonly out_dir="data/simulation/mt"
readonly out_prefix="${out_dir}/ukb_eur_25k_samples_chr${chr}"

readonly hap_file="/well/lindgren/flassen/ressources/hapmap/weights.l2.ldscore.liftover.ht"
readonly annotation_table="data/vep/hail/ukb_wes_200k_chr${chr}_vep.ht"
readonly sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.bed" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --in_prefix "${in_prefix}"\
     --in_type "mt" \
     --random_samples 25000 \
     --sample_seed 42 \
     --out_prefix "${out_prefix}" \
     --out_type "mt"
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi




