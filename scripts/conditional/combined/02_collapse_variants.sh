#!/usr/bin/env bash
#
# @description Append pseudo variants with actual variants for downstream conditional analysis.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=collapse_variants
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/collapse_variants.log
#SBATCH --error=logs/collapse_variants.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 4
#SBATCH --array=1-22


set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/combined/02_collapse_variants.R"

readonly chr="${SLURM_ARRAY_TASK_ID}"

readonly vep_dir="data/mt/prefilter/final_90_loftee"
readonly vep="${vep_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.csqs.txt.gz"

readonly rare_dir="data/conditional/rare/combined/mt"
readonly rare_file="${rare_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_markers.txt.gz"

readonly vcf_dir="data/conditional/rare/combined/mt"
readonly vcf="${vcf_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense.vcf.bgz"

readonly out_dir="data/conditional/combined/collapsed"
readonly out_prefix="ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_collapsed/sleep"
readonly out="${out_prefix}.vcf"

mkdir -p ${out_dir}

if [ ! -f "${out}" ]; then
  set_up_rpy
  Rscript "${rscript}" \
    --path_worst_csq_by_gene_canonical ${vep} \
    --path_markers ${rare_file} \
    --path_vcf ${vcf} \
    --chrom ${chr} \
    --out_prefix ${out_prefix}
fi

if [ ! -f "${out}.gz" ]; then
  module load BCFtools/1.12-GCC-10.3.0
  bgzip -c ${out} > "${out}.gz"
fi

if [ ! -f "${out}.gz.csi" ]; then
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out}.gz" "csi"
fi





