#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_alt_alleles
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/export_alt_alleles.log
#SBATCH --error=logs/export_alt_alleles.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-21

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/vcf_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

#readonly in_dir="data/mt/prefilter/pp90"
#readonly in_vcf="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.vcf.gz"

#readonly out_dir="data/mt/prefilter/pp90"
#readonly out_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.alt_alleles"

readonly in_dir="data/mt/prefilter/no_pp_cutoff/old"
readonly in_vcf="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.vcf.bgz"

readonly out_dir="data/mt/prefilter/no_pp_cutoff/old"
readonly out_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.alt_alleles"




module load BCFtools
bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT %PP %AC %AN %AF\n]' ${in_vcf}\
  |  awk '{print $1" "$2" "$3" "$4" "$5" "$5" "$6}' | gzip > ${out_prefix}.txt.gz




