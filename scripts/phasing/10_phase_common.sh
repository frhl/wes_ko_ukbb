#!/usr/bin/env bash
#
# @description phase genotyping array calls using SHAPEIT5
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=phase_common
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/phase_common.log
#SBATCH --error=logs/phase_common.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 16
#SBATCH --array=20-22
#
#$ -N phase_common
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phase_common.log
#$ -e logs/phase_common.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 22
#$ -V

# 16 cores seems to be sufficient ~8-12H

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly in_dir="data/unphased/calls/prefilter/500k"
readonly out_dir="data/phased/calls/shapeit5/500k"
readonly ref_dir="/well/lindgren/flassen/ressources/panels/liftover_reference_panel/data/liftover"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_file="${in_dir}/ukb_split_calls_500k_chr${chr}_no_parents.vcf.bgz"
readonly out_prefix="${out_dir}/ukb_phased_calls_500k_chr${chr}"
readonly out="${out_prefix}.vcf.gz"

readonly ref="${ref_dir}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.bgz"
readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"
readonly threads=$( get_threads )

mkdir -p ${out_dir}

if [ ! -f ${out} ]; then
  SECONDS=0
  set_up_shapeit5
  ${SHAPEIT_phase_common} \
    --input ${in_file} \
    --map ${gmap} \
    --region "chr${chr}" \
    --thread ${threads} \
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





