#!/usr/bin/env bash
#
# @description phase genotyping array calls using SHAPEIT5
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=phase_rare
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/phase_rare.log
#SBATCH --error=logs/phase_rare.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=21
#
#$ -N phase_rare
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phase_rare.log
#$ -e logs/phase_rare.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly ref_dir="/well/lindgren/flassen/ressources/panels/liftover_reference_panel/data/liftover"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly common_dir="data/phased/calls/shapeit5/200k_from_500k"
readonly common_vcf="${common_dir}/ukb_phased_calls_200k_from_500k_chr${chr}.vcf.bgz"

readonly rare_dir="data/unphased/wes_union_calls/200k/test_with_phased_calls"
#readonly rare_dir="data/unphased/wes_union_calls/bcftools/newtest"
#readonly rare_vcf="${rare_dir}/ukb_split_wes_200k_chr${chr}_no_parents.vcf.bgz"
readonly rare_vcf="${rare_dir}/ukb_wes_union_calls_chr${chr}.vcf.gz"

readonly out_dir="data/phased/wes_scaffold_calls/test18"
readonly out_prefix="${out_dir}/ukb_shapeit5_full_200k_from_500k_chr${chr}"
readonly out="${out_prefix}.vcf.gz"
readonly log="${out_prefix}.log"

readonly ref="${ref_dir}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.bgz"
readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"
readonly threads=$( get_threads )

readonly population="190000"

# Note, we need to double check the sample ordering of these files. 
vcf_check_sample_order ${rare_vcf} ${common_vcf}
echo "CALLS: ${common_vcf}"
echo WES: ${rare_vcf}""

#readonly region="chr${chr}:32413059-32437972"
readonly region="chr${chr}:12413059-42437972"

mkdir -p ${out_dir}

if [ ! -f ${out} ]; then
  SECONDS=0
  set_up_shapeit5
  ${SHAPEIT_phase_rare} \
    --input-plain ${rare_vcf} \
    --input-region ${region} \
    --scaffold ${common_vcf} \
    --scaffold-region ${region} \
    --thread ${threads} \
    --output ${out} \
    --log ${log} \
    && print_update "Finished phasing variants for chr${chr}, out: ${out}" "${SECONDS}" \
    || raise_error "$( print_update "Phasing variants failed for chr${chr}" ${SECONDS} )"
    module purge
else
  >&2 echo "${out} already exists. Skipping"  
fi

if [ ! -f "${out}.tbi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out}" "tbi"
fi





