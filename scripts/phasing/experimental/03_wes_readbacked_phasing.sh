#!/usr/bin/env bash
#
# @description phase variants from UKBB whole exome sequencing using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_readbacked_phasing
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_readbacked_phasing.log
#SBATCH --error=logs/wes_readbacked_phasing.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=21


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/unphased/wes/prefilter"
readonly out_dir="data/phased/wes/naive"

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly in_file="${in_dir}/ukb_eur_wes_prefilter_200k_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_eur_wes_prefilter_200k_naive_phasing_chr${chr}.vcf.gz"
readonly fam_file="${fam_dir}/ukb11867_pedigree.fam"

# human genome reference
readonly grch38="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly cram_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly bam_dir="data/reads/bam"
# placeholder files to be subsetted with sample ID
readonly cram_placeholder="${cram_dir}/SAMPLE_oqfe.cram"
readonly bam_placeholder="${bam_dir}/SAMPLE_oqfe.bam"


cram_to_bam() {
  # load samtools
  module load SAMtools/1.13-GCC-10.3.0
  # convert cram to bam file  
  local _reference=${1}
  local _cram=${2}
  local _bam=${3}
  samtools view -b -T ${_reference} -o ${_bam} ${_cram}
}

readonly eid="1000278"
readonly cram_file=${cram_placeholder/SAMPLE/${eid}}
readonly bam_file=${bam_placeholder/SAMPLE/${eid}}
cram_to_bam ${grch38} ${cram_file} ${bam_file}












