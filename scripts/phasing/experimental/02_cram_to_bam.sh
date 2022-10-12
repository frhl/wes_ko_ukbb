#!/usr/bin/env bash
#
# @description phase variants from UKBB whole exome sequencing using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=cram_to_bam
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/cram_to_bam.log
#SBATCH --error=logs/cram_to_bam.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1

# require to run the following
module load samtools/1.8-gcc5.4.0

# human genome reference
readonly grch38="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly cram_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly bam_dir="data/reads/bam"
# placeholder files to be subsetted with sample ID
readonly cram_placeholder="${cram_dir}/SAMPLE_oqfe.cram"
readonly bam_placeholder="${bam_dir}/SAMPLE_oqfe.bam"

cram_to_bam() {
  # load samtools
  # convert cram to bam file  
  local _reference=${1}
  local _cram=${2}
  local _bam=${3}
  if [ ! -f ${_bam} ]; then
    samtools view -b -T ${_reference} -o ${_bam} ${_cram}
  else
    echo "${bam_file} already exists. Skipping.."
  fi
  samtools index ${_bam}
}

readonly eid="2183281" #"5101274" #"1000278"
readonly cram_file=${cram_placeholder/SAMPLE/${eid}}
readonly bam_file=${bam_placeholder/SAMPLE/${eid}}
cram_to_bam ${grch38} ${cram_file} ${bam_file}

