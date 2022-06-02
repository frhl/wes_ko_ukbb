#!/usr/bin/env bash
#
#$ -N calls_naive_phasing
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/calls_naive_phasing.log
#$ -e logs/calls_naive_phasing.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 6
#$ -q short.qc@@short.hga
#$ -t 20-22


source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly in_dir="data/unphased/calls"
readonly out_dir="data/phased/calls/shapeit4/naive"
readonly ref_dir="/well/lindgren/flassen/ressources/panels/liftover_reference_panel/data/liftover"

readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"

readonly chr="${SGE_TASK_ID}"
readonly in_file="${in_dir}/ukb_prefilter_calls_200k_chr${chr}.vcf.bgz"
readonly out_prefix="${out_dir}/ukb_prefilter_calls_200k_chr${chr}"
readonly out="${out_prefix}.vcf.gz"
readonly trio="${out_prefix}.trio"

readonly ref="${ref_dir}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.bgz"
readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"

mkdir -p ${out}

if [ ! -f ${out} ]; then
  module load SHAPEIT4/4.2.2-foss-2021a
  SECONDS=0
  shapeit4.2 \
    --input ${in_file} \
    --map ${gmap} \
    --region "chr${chr}" \
    --thread $(( ${NSLOTS}-1 )) \
    --output ${out} \
    && print_update "Finished phasing variants for chr${chr}, out: ${out}" "${SECONDS}" \
    || raise_error "$( print_update "Phasing variants failed for chr${chr}" ${SECONDS} )"
    module purge
fi

if [ ! -f "${out}.tbi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out}" "tbi"
fi

if [ ! -f ${trio} ]; then
    bcftools +trio-switch-rate ${out} -- -p ${pedigree} > ${trio}
    switch_errors_by_site ${out} ${pedigree}
fi





