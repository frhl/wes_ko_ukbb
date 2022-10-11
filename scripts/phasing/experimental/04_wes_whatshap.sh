#!/usr/bin/env bash
#
# @description phase variants from UKBB whole exome sequencing using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_whatshap
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_whatshap.log
#SBATCH --error=logs/wes_whatshap.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

#readonly in_dir="data/mt/annotated/"
readonly in_dir="data/unphased/wes/prefilter"
readonly out_dir="data/phased/experimental"

mkdir -p ${out_dir}

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly vcf_file="${in_dir}/ukb_eur_wes_prefilter_200k_chr21.vcf.bgz"
#readonly vcf_file="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_eur_wes_200k_chr${chr}_whatshap_phased"

# human genome reference
readonly grch38="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly cram_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly cram_placeholder="${cram_dir}/SAMPLE_oqfe.cram"

readonly eid="1281289" #"1000278"
readonly cram_file=${cram_placeholder/SAMPLE/${eid}}
readonly bam_file=${bam_placeholder/SAMPLE/${eid}}

set_up_whatshap
whatshap phase \
  --reference ${grch38} \
  --chromosome chr${chr} \
  --sample  ${eid} \
  -o ${out_file} \
  ${vcf_file} \
  ${cram_file}






