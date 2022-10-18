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
#SBATCH --array=1


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

#readonly in_dir="data/mt/annotated/"
readonly in_dir="data/unphased/wes/prefilter"
readonly out_dir="data/phased/experimental"

mkdir -p ${out_dir}

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly vcf_file="${in_dir}/ukb_eur_wes_prefilter_200k_chr${chr}.vcf.bgz"
#readonly vcf_file="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}_phased"

# human genome reference
readonly grch38="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly cram_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly bam_dir="data/reads/bam"
# placeholder files to be subsetted with sample ID
readonly cram_placeholder="${cram_dir}/SAMPLE_oqfe.cram"
readonly bam_placeholder="${bam_dir}/SAMPLE_oqfe.bam"

# SAMtools is also required for PHASER to run
#module load SAMtools/1.13-GCC-10.3.0
module load samtools/1.8-gcc5.4.0
module load htslib/1.8-gcc5.4.0 
module load BEDTools/2.29.2-GCC-8.3.0
module load BCFtools/1.10.2-GCC-8.3.0

cram_to_bam() {
  # load samtools
  # convert cram to bam file  
  local _reference=${1}
  local _cram=${2}
  local _bam=${3}
  samtools view -b -T ${_reference} -o ${_bam} ${_cram}
  samtools index ${_bam}
}

readonly eid="5101274" #"1000278"
readonly cram_file=${cram_placeholder/SAMPLE/${eid}}
readonly bam_file=${bam_placeholder/SAMPLE/${eid}}
if [ ! -f ${bam_file} ]; then
  cram_to_bam ${grch38} ${cram_file} ${bam_file}
else
  echo "${bam_file} already exists. Skipping.."
fi

readonly threads=42

# use read backed phasing
#set -x
#set_up_phaser
#python2.7 ${PHASER_PATH} \
#   --vcf ${vcf_file} \
#   --bam ${bam_file} \
#   --threads ${threads} \
#   --pass_only 0 \
#   --chr "chr${chr}" \
#   --baseq 30 \
#   --mapq 15 \
#   --paired_end 1 \
#   --sample ${eid} \
#   --o ${out_file}
#set +x







