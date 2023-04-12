#!/usr/bin/env bash
#
# Append pseudo variants with actual variants for downstream conditional analysis.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=combine_additive_recessive
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/combine_additive_recessive.log
#SBATCH --error=logs/combine_additive_recessive.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/dominance/01_combine_additive_recessive.py"

readonly array_idx=$( get_array_task_id )
readonly chr="$( get_chr ${array_idx} )"

# note, that common/rare variants are already included in recessive encoding
readonly additive_dir="data/knockouts/alt/pp90/encoding_012"
readonly recessive_dir="data/conditional/combined/combine_collapsed_urv"

# additive encoding 0 1 2
readonly additive_path_wo_ext="${additive_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense"
readonly additive_path="${additive_path_wo_ext}.vcf.bgz"
readonly additive_type="vcf"

# recessive/dominance encoding 0 0 1
readonly recessive_path_wo_ext="${recessive_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense"
readonly recessive_path="${recessive_path_wo_ext}.vcf.bgz"
readonly recessive_type="vcf"

readonly out_dir="data/conditional/dominance/combine_encodings"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense"
readonly out_type="vcf"
readonly out="${out_prefix}.vcf.bgz"

mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --chrom ${chr} \
     --additive_path ${additive_path} \
     --additive_type ${additive_type} \
     --recessive_path ${recessive_path} \
     --recessive_type ${recessive_type} \
     --out_prefix ${out_prefix} \
     --out_type ${out_type}
fi

# index the resulting file.
if [ ! -f "${out_prefix}.vcf.csi" ] & [ "${out_type}" == "vcf" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix}.vcf.bgz" "csi"
fi



