#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=annotate_ps
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/annotate_ps.log
#SBATCH --error=logs/annotate_ps.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 3
#SBATCH --array=21
#
#$ -N annotate_ps
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/annotate_ps.log
#$ -e logs/annotate_ps.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly hail_script="scripts/phasing/phasing/06_annotate_ps.py"
readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly ref_dir="data/prephased/wes_union_calls"
readonly ref_path="${ref_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.gz"
readonly ref_type="vcf"

readonly phased_dir="data/phased/wes_union_calls/200k/shapeit5/phase_rare/ukb_wes_union_calls_shapeit5_200k_chr${chr}-20xshort"
readonly phased_path="${phased_dir}/shapeit5_prs100000_pro25000_mprs150000.1of1.vcf.gz"
readonly phased_type="vcf"

readonly out_dir="data/phased/wes_union_calls/200k/calibration"
readonly out_prefix="${out_dir}/ukb_shapeit5_whatshap_chr${chr}"
readonly out_type="mt"

mkdir -p ${out_dir}


# combine parents and children in same vcf
# note: we can't combine data using bcftools beacuse some variants
# share the same position (but not same reference alleles)
module purge
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --phased_path ${phased_path} \
  --phased_type ${phased_type} \
  --ref_path ${ref_path} \
  --ref_type ${ref_type} \
  --out_prefix ${out_prefix} \
  --out_type ${out_type}

