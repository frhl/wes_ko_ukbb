#!/usr/bin/env bash
#
# @description run phasing on unphased singletons using informative variants
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=readbacked_phasing
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/readbacked_phasing.log
#SBATCH --error=logs/readbacked_phasing.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-22


source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/experimental/03_readbacked_phasing.py"
readonly reference="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
# SNPs that are to be assed
readonly query_snps_dir="data/reads/singletons"
readonly query_snps="${query_snps_dir}/samples_with_damaging_missense_singletons.txt"
# Nearby SNPs used as reference points for phasing
readonly informative_snps_dir="data/mt/heterozygotes/improved" 
readonly informative_snps="${informative_snps_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.tsv.gz"
# outfiles
readonly out_dir="data/reads/phased"
readonly outfile="${out_dir}/ukb_eur_wes_singletons_chr${chr}.txt"

mkdir -p ${out_dir}

# run phasing
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --path_fasta_reference ${reference} \
   --path_informative_snps ${informative_snps} \
   --path_query_snps ${query_snps} \
   --chromosome "chr${chr}" \
   --outfile ${outfile}




