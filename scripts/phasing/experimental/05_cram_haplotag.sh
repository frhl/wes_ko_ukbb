#!/usr/bin/env bash
#
# @description phase variants from UKBB whole exome sequencing using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=cram_haplotag
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/cram_haplotag.log
#SBATCH --error=logs/cram_haplotag.errors.log
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
readonly vcf_file="${in_dir}/eid1281289_eur_wes_union_calls_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_eur_wes_200k_chr${chr}_haplotag.vcf"

#readonly grch38="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly cram_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly bam_dir="data/reads/bam"
readonly cram_placeholder="${cram_dir}/SAMPLE_oqfe.cram"
readonly bam_placeholder="${bam_dir}/SAMPLE_oqfe.bam"

readonly eid="1281289"
readonly cram_file=${cram_placeholder/SAMPLE/${eid}}
readonly bam_file=${bam_placeholder/SAMPLE/${eid}}
readonly threads=30

echo "# REFERENCE: ${grch38}"
echo "# BAM FILE: ${bam_file}"
echo "# CRAM FILE: ${cram_file}"
echo "# VCF FILE: ${vcf_file}"

module load htslib/1.8-gcc5.4.0


readonly refdir="/well/lindgren/flassen/ressources/genome_reference/broad"
readonly grch38="${refdir}/Homo_sapiens_assembly38.fasta"
export REF_PATH="${refdir}/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="${refdir}/ref/cache/%2s/%2s/%s"

SECONDS=0
set_up_whatshap
set -x
whatshap haplotag \
  --reference="${grch38}" \
  --output="${out_file}" \
  --sample="${eid}" \
  --ignore-read-groups \
  ${vcf_file} \
  ${cram_file}
set +x
echo "Time elapsed ${SECONDS}"





