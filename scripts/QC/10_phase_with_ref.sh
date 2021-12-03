#!/usr/bin/env bash
#
#$ -N refphase
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/refphase.log
#$ -e logs/refphase.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qc@@short.hge
#$ -t 21

source utils/qsub_utils.sh

readonly in_dir="data/unphased/post-qc/"
readonly out_dir="data/phased/test-phasing"

readonly chr="${SGE_TASK_ID}"
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_wes_200k_refphased_chr${chr}.vcf.gz"

readonly ref=""
readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"

mkdir -p ${out_dir}
module load SHAPEIT4/4.2.0-foss-2019b
shapeit4 --input ${OUT_ROOT}${IN_VCF} \
  --map ${gmap} \
  --region ${chr} \
  --thread 20 \
  --reference ${ref} \
  --output ${out_file}



