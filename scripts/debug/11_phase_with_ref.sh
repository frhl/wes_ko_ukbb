#!/usr/bin/env bash
#
#$ -N refphase
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/refphase.log
#$ -e logs/refphase.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 24
#$ -q long.qc@@long.hga
#$ -t 21

source utils/qsub_utils.sh
source utils/bash_utils.sh

#readonly in_dir="data/unphased/post-qc"
readonly in_dir="data/phased/test-phasing"
readonly out_dir="data/phased/test-phasing"
readonly ref_dir="/well/lindgren/flassen/ressources/panels/liftover_reference_panel/data/liftover"

readonly chr="${SGE_TASK_ID}"
#readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
readonly in_file="${in_dir}/ukb_wes_200k_filtered_maf0001_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_wes_200k_refphased_chr${chr}.vcf.gz"

readonly ref="${ref_dir}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.bgz"
readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"

mkdir -p ${out_dir}
module load SHAPEIT4/4.2.2-foss-2021a
SECONDS=0
shapeit4.2 \
  --input ${in_file} \
  --map ${gmap} \
  --region "chr${chr}" \
  --thread $(( ${NSLOTS}-1 )) \
  --reference ${ref} \
  --output ${out_file} \
  --sequencing \
  && print_update "Finished phasing variants for chr${chr}, out: ${out}" "${SECONDS}" \
  || raise_error "$( print_update "Phasing common variants failed for chr${chr}" ${SECONDS} )"



