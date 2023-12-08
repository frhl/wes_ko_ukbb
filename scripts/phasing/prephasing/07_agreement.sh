#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=agreement
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/agreement.log
#SBATCH --error=logs/agreement.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 5

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly rscript="scripts/phasing/prephasing/07_agreement.R"

readonly input_dir="data/prephased/wes_union_calls/10k" # 10k samples all autosomes
readonly input_path="${input_dir}/ukb_shapeit5_whatshap_chrCHR.10k.PP.PS.txt.gz"

readonly out_dir="data/prephased/wes_union_calls/10k"
readonly out_prefix="${out_dir}/ukb_shapeit5_whatshap.10k.PP.PS.agreement"

#readonly input_dir="data/prephased/wes_union_calls/full_phase_conf" # all chromosomes 1000 samples
#readonly input_path="${input_dir}/ukb_shapeit5_whatshap_chrCHR.PP.PS.txt.gz"

#readonly out_dir="data/prephased/wes_union_calls/full_phase_conf"
#readonly out_prefix="${out_dir}/ukb_shapeit5_whatshap.PP.PS.agreement"

readonly wes_variants="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"
readonly ko_samples="data/phenotypes/samples/ukb_wes_ko.qc.nfe.samples"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --input_path "${input_path}" \
  --sites "${wes_variants}" \
  --samples "${ko_samples}" \
  --out_prefix "${out_prefix}"



