#!/usr/bin/env bash
#
# @description phase variants from UKBB whole exome sequencing using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=cram_multiphase
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/cram_multiphase.log
#SBATCH --error=logs/cram_multiphase.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh


#readonly in_dir="data/reads/samples"


readonly chr="${SLURM_ARRAY_TASK_ID}"

readonly in_dir="data/unphased/wes_union_calls"
readonly vcf_file="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"

readonly out_dir="data/reads/phased"
readonly out_file="${out_dir}/eid1281289_5101274_eur_wes_union_calls_chr${chr}_phased_full_vcf.vcf"

#readonly vcf_file="${in_dir}/eid1281289_5101274_eur_wes_union_calls_chr${chr}.vcf.bgz"
#readonly out_file="${out_dir}/eid1281289_5101274_eur_wes_union_calls_chr${chr}_phased.vcf"



#readonly grch38="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly cram_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly bam_dir="data/reads/bam"
readonly cram_placeholder="${cram_dir}/SAMPLE_oqfe.cram"

readonly cram_file1=${cram_placeholder/SAMPLE/1281289}
readonly cram_file2=${cram_placeholder/SAMPLE/5101274}
module load htslib/1.8-gcc5.4.0

readonly refdir="/well/lindgren/flassen/ressources/genome_reference/broad"
readonly grch38="${refdir}/Homo_sapiens_assembly38.fasta"
export REF_PATH="${refdir}/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="${refdir}/ref/cache/%2s/%2s/%s"

mkdir -p ${out_dir}

SECONDS=0
set_up_whatshap
whatshap phase \
  --reference="${grch38}" \
  --output="${out_file}" \
  --indels \
  ${vcf_file} \
  ${cram_file1} \
  ${cram_file2}
echo "Time elapsed ${SECONDS}"





