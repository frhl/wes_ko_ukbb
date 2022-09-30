#!/usr/bin/env bash
#
# @description filter WES quality-controlled data.
# @note singletons are excluded here but re-introduced in a later script.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prefilter_wes
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prefilter_wes.log
#SBATCH --error=logs/prefilter_wes.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=20-22


source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/01_prefilter_wes.py"
readonly in_dir="data/unphased/wes/post-qc"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/unphased/wes/prefilter"
readonly out_prefix="${out_dir}/ukb_eur_wes_prefilter_200k_chr${chr}"
readonly out_type="vcf"

mkdir -p ${out_dir}
if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --input_path "${in_file}" \
     --input_type "${in_type}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --exclude_trio_parents \
     --ancestry "eur" \
     --min_mac 2 \
     --missing 0.05
  set +x
  log_runtime ${SECONDS}
else
  print_update "file ${out} already exists. Skipping!"
fi

module purge
module load BCFtools/1.12-GCC-10.3.0
make_tabix "${out_prefix}.vcf.bgz" "tbi"



