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


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/unphased/wes/prefilter"
readonly out_dir="data/phased/experimental"

mkdir -p ${out_dir}

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly vcf_file="${in_dir}/ukb_eur_wes_prefilter_200k_chr21.vcf.bgz"
#readonly vcf_file="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_eur_wes_200k_chr${chr}_smartphase.tsv"

# all variants in probabilistic knockouts
# allele counts of variants

# variants to be assessed
readonly variants="data/reads/variants/test"
# human genome reference
readonly grch38="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly cram_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly bam_dir="data/reads/bam"
# placeholder files to be subsetted with sample ID
readonly cram_placeholder="${cram_dir}/SAMPLE_oqfe.cram"
readonly bam_placeholder="${bam_dir}/SAMPLE_oqfe.bam"
# select specific sample
readonly eid="5101274"
readonly cram_file=${cram_placeholder/SAMPLE/${eid}}
readonly bam_file=${bam_placeholder/SAMPLE/${eid}}
readonly mapq=30

set_up_smartphase

set -x
SECONDS=0
java -jar ${SMARTPHASE_PATH} \
  -f ${variants} \
  -p ${eid} \
  -r ${bam_file} \
  -m ${mapq} \
  -a ${vcf_file} \
  -o ${out_file}
echo $SECONDS
set +x




