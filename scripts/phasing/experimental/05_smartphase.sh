#!/usr/bin/env bash
#
# @description phase variants from UKBB whole exome sequencing using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=smartphase
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/smartphase.log
#SBATCH --error=logs/smartphase.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=21

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/mt/annotated"
readonly var_dir="data/reads/singletons"
readonly out_dir="data/reads/phased"

mkdir -p ${out_dir}

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly vcf_file="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}_smartphase.tsv"
readonly var_file="${var_dir}/samples_with_damaging_missense_singletons.txt"


# all variants in probabilistic knockouts
# allele counts of variants

# human genome reference
readonly grch38="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly cram_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly bam_dir="data/reads/bam"
# placeholder files to be subsetted with sample ID
readonly cram_placeholder="${cram_dir}/SAMPLE_oqfe.cram"
readonly bam_placeholder="${bam_dir}/SAMPLE_oqfe.bam"
# select specific sample
readonly eid=2613995 #"5101274"
readonly cram_file=${cram_placeholder/SAMPLE/${eid}}
readonly bam_file=${bam_placeholder/SAMPLE/${eid}}
readonly mapq=30

# get singletons save at a specific file
readonly variants="${out_dir}/input_${eid}.csv"
cat ${var_file} | grep -E "chr${chr}:" | grep "${eid}" | cut -f4 > "${variants}"

if [ $( cat ${variants} | wc -l ) -gt 0 ]; then
  if [ -f ${bam_file} ]; then
  set_up_smartphase
  set -x
  java -jar ${SMARTPHASE_PATH} \
    --all-variants ${vcf_file} \
    --filtered-variants ${variants} \
    --patient ${eid} \
    --reads ${bam_file} \
    --mapq ${mapq} \
    --output ${out_file} \
    --cutoff 0.1 \
    --reject-phase \
    -vcf
  set +x
  else
    >&2 echo "${bam_file} does not exist!"
  fi
else 
  >&2 echo "No variants for sample ${eid}"
fi



