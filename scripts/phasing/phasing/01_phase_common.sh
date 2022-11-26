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
#SBATCH --cpus-per-task 18
#SBATCH --array=3,5,9
#
#$ -N phase_common
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phase_common.log
#$ -e logs/phase_common.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -q short.qc
#$ -t 21
#$ -V

# For 500K, chr21 takes ~23 Hours.
# Will need to use the long queue for the rest

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly tranche="200k"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/unphased/wes_union_calls/prefilter_no_maf_cutoff/${tranche}"
readonly in_file="${in_dir}/ukb_wes_union_calls_chr${chr}.vcf.gz"

readonly out_dir="data/phased/wes_union_calls/${tranche}/shapeit5/phase_common/newrun"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_${tranche}_chr${chr}_phase_common"
readonly out="${out_prefix}.vcf.gz"
readonly log="${out_prefix}.log"

readonly ref_dir="/well/lindgren/flassen/ressources/panels/liftover_reference_panel/data/liftover"
readonly ref="${ref_dir}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.bgz"
readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"
readonly threads=$( get_threads )

readonly pbwt_modulo=0.1 # deault is 0.1, in shapeit4 with sequencing it is 0.0002
readonly pbwt_depth=4 # deault is 4
readonly pbwt_mac=5 # deafult is 5
readonly pbwt_mdr=0.1 # default is 0.1
readonly min_maf=0.001 # no default value

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
    --pbwt-modulo ${pbwt_modulo} \
    --pbwt-depth ${pbwt_depth} \
    --pbwt-mac ${pbwt_mac} \
    --pbwt-mdr ${pbwt_mdr} \
    --filter-maf ${min_maf} \
    --log "${log}" \
    && print_update "Finished phasing variants for chr${chr}, out: ${out}" "${SECONDS}" \
    || raise_error "$( print_update "Phasing variants failed for chr${chr}" ${SECONDS} )"
    module purge
fi

if [ ! -f "${out}.tbi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out}" "tbi"
fi





