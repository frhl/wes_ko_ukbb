#!/usr/bin/env bash
#
# @description phase combined set of UK Biobank whole exome and genotyping array using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_union_calls_phasing
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_union_calls_phasing.log
#SBATCH --error=logs/wes_union_calls_phasing.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --array=21
#
#$ -N wes_union_calls_phasing
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_union_calls_phasing.log
#$ -e logs/wes_union_calls_phasing.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 16
#$ -q short.qc
#$ -t 21
#$ -V

# note, shmem 20 is not enough for chrom 21.

source utils/qsub_utils.sh
source utils/vcf_utils.sh
source utils/bash_utils.sh


readonly in_dir="data/unphased/wes_union_calls"
#readonly in_dir="data/prephased/wes_union_calls"
readonly out_dir="data/phased/wes_union_calls/benchmark/200k"
readonly ref_dir="/well/lindgren/flassen/ressources/panels/liftover_reference_panel/data/liftover"
readonly fam_dir="/well/lindgren/UKBIOBANK/nbaya/resources"

readonly chr="$( get_array_task_id )"
readonly in_file="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.gz"
readonly out_file="${out_dir}/ukb_wes_union_calls_200k_phased_shapeit4_chr${chr}.vcf.gz"

readonly ref="${ref_dir}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.bgz"
readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"

mkdir -p ${out_dir}

if [ ! -f ${out_file} ]; then
  module load SHAPEIT4/4.2.2-foss-2021a
  SECONDS=0
  shapeit4.2 \
    --input ${in_file} \
    --map ${gmap} \
    --region "chr${chr}" \
    --thread 15 \
    --pbwt-mac 2 \
    --output ${out_file} \
    --use-PS 0.0001 \
    --sequencing \
    && print_update "Finished phasing variants for chr${chr}, out: ${out_file}" "${SECONDS}" \
    || raise_error "$( print_update "Phasing variants failed for chr${chr}" ${SECONDS} )"
  module purge
fi

if [ ! -f "${out_file}.tbi" ]; then
     module load BCFtools/1.12-GCC-10.3.0
     make_tabix "${out_file}" "tbi"
fi



