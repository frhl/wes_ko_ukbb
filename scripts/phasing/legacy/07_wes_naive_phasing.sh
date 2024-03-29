#!/usr/bin/env bash
#
# @description phase variants from UKBB whole exome sequencing using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_naive_phasing
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_naive_phasing.log
#SBATCH --error=logs/wes_naive_phasing.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 10
#SBATCH --array=21


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/unphased/wes/from_wes_union_calls"
readonly out_dir="data/phased/wes/from_wes_union_calls"
readonly ref_dir="/well/lindgren/flassen/ressources/panels/liftover_reference_panel/data/liftover"

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly in_file="${in_dir}/ukb_wes_200k_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_wes_200k_naive_phased_chr${chr}.vcf.gz"

readonly ref="${ref_dir}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.bgz"
readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"
readonly threads=$(( ${SLURM_CPUS_ON_NODE} - 1))

mkdir -p ${out_dir}

if [ ! -f "${out_file}" ]; then
  module load SHAPEIT4/4.2.2-foss-2021a
  SECONDS=0
  shapeit4.2 \
    --input ${in_file} \
    --map ${gmap} \
    --region "chr${chr}" \
    --thread ${threads} \
    --output ${out_file} \
    --sequencing \
    && print_update "Finished phasing variants for chr${chr}, out: ${out}" "${SECONDS}" \
    || raise_error "$( print_update "Phasing variants failed for chr${chr}" ${SECONDS} )"
fi


if [ ! -f "${out_file}.tbi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_file}" "tbi"
fi

