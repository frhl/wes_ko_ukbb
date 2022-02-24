#!/usr/bin/env bash
#
#$ -N wes_naive_phasing
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_naive_phasing.log
#$ -e logs/wes_naive_phasing.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 20

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/unphased/wes/post-qc"
readonly out_dir="data/phased/wes/naive/"
readonly ref_dir="/well/lindgren/flassen/ressources/panels/liftover_reference_panel/data/liftover"
readonly fam_dir="/well/lindgren/UKBIOBANK/nbaya/resources"

readonly chr="${SGE_TASK_ID}"
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_wes_200k_naive_phasing_chr${chr}.vcf.gz"
readonly ser_file="${out_dir}/ukb_wes_200k_naive_phasing_chr${chr}.txt"
readonly fam_file="${fam_dir}/ukb11867_pedigree.fam"

readonly ref="${ref_dir}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.bgz"
readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"

mkdir -p ${out_dir}

if [ ! -f "${out_file}" ]; then
  module load SHAPEIT4/4.2.2-foss-2021a
  SECONDS=0
  shapeit4.2 \
    --input ${in_file} \
    --map ${gmap} \
    --region "chr${chr}" \
    --thread $(( ${NSLOTS}-1 )) \
    --output ${out_file} \
    --sequencing \
    && print_update "Finished phasing variants for chr${chr}, out: ${out}" "${SECONDS}" \
    || raise_error "$( print_update "Phasing variants failed for chr${chr}" ${SECONDS} )"
  module purge
  log_runtime ${SECONDS}
fi

module load BCFtools/1.12-GCC-10.3.0

if [ ! -f "${out_file}.tbi" ]; then
  module purge
  make_tabix "${out_file}" "tbi"
fi

export BCFTOOLS_PLUGINS="/apps/eb/2020b/skylake/software/BCFtools/1.12-GCC-10.3.0/libexec/bcftools"
bcftools +trio-switch-rate "${out_file}" -- -p "${fam_file}" > "${ser_file}"


