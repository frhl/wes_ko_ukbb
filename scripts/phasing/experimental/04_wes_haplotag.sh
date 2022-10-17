#!/usr/bin/env bash
#
# @description phase variants from UKBB whole exome sequencing using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_haplotag
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_haplotag.log
#SBATCH --error=logs/wes_haplotog.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh


#readonly in_dir="data/mt/annotated/"
readonly in_dir="data/reads/samples"
readonly out_dir="data/reads/haplotag"

mkdir -p ${out_dir}

readonly chr="${SLURM_ARRAY_TASK_ID}"
#readonly vcf_file="${in_dir}/ukb_eur_wes_prefilter_200k_chr21.vcf.bgz"
readonly vcf_file="${in_dir}/ukb_wes_union_calls_200k_chr21.vcf.bgz"
readonly out_file="${out_dir}/ukb_eur_wes_200k_chr${chr}_whatshap_phased"

readonly grch38="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly cram_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly cram_placeholder="${cram_dir}/SAMPLE_oqfe.cram"

readonly eid="1281289" #"1000278"
readonly cram_file=${cram_placeholder/SAMPLE/${eid}}
readonly threads=30

SECONDS=0
set_up_whatshap
whatshap haplotag \
  --reference ${grch38} \
  --chromosome chr${chr} \
  --output-threads ${threads} \
  --sample  ${eid} \
  -o ${out_file} \
  ${vcf_file} \
  ${cram_file}
echo "Time elapsed ${SECONDS}"





