#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_ps_pp
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/export_ps_pp.log
#SBATCH --error=logs/export_ps_pp.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly hail_script="scripts/phasing/prephasing/06_export_ps_pp.py"
readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/prephased/wes_union_calls/revision"
readonly in_path="${in_dir}/ukb_shapeit5_whatshap_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/prephased/wes_union_calls/50k"
readonly out_prefix="${out_dir}/ukb_shapeit5_whatshap_chr${chr}.50k"

mkdir -p ${out_dir}

# read-backed phasing file with the samples that have been processed
readonly ref_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/readbacked/data/phased/subset/50k"
readonly ref_path="${ref_dir}/UKB.wes.200k.chr${chr}.50k.b1of4.vcf.gz"
readonly ref_type="vcf"

# samples to subset from that has read-backed phasing
readonly sample_file="${out_prefix}.samples"
module load BCFtools
#bcftools query -l ${ref_path} | head -n 100000 > ${sample_file}
bcftools query -l ${ref_path} > ${sample_file}


module purge
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --sample_file ${sample_file} \
  --phased_path ${in_path} \
  --phased_type ${in_type} \
  --out_prefix ${out_prefix}

rm ${sample_file}

