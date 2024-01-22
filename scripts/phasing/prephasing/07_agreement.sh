#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=agreement
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/agreement.log
#SBATCH --error=logs/agreement.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 2

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly rscript="scripts/phasing/prephasing/07_agreement.R"

readonly input_dir="data/prephased/wes_union_calls/50k" # 10k samples all autosomes
readonly input_path="${input_dir}/ukb_shapeit5_whatshap_chrCHR.50k.PP.PS.txt.gz"

readonly out_dir="data/prephased/wes_union_calls/50k"
readonly out_prefix="${out_dir}/ukb_shapeit5_whatshap.50k.PP.PS.agreement.long"

#readonly input_dir="data/prephased/wes_union_calls/full_phase_conf" # all chromosomes 1000 samples
#readonly input_path="${input_dir}/ukb_shapeit5_whatshap_chrCHR.PP.PS.txt.gz"

#readonly out_dir="data/prephased/wes_union_calls/full_phase_conf"
#readonly out_prefix="${out_dir}/ukb_shapeit5_whatshap.PP.PS.agreement"

readonly wes_variants="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

readonly samples_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/samples"

#readonly ko_samples="data/phenotypes/samples/ukb_wes_ko.qc.nfe.samples"

mkdir -p ${out_dir}
set_up_rpy

for anc in "eur"; do
  samples="${samples_dir}/UKB.chr21.samples.${anc}.txt"
  out_prefix_anc="${out_prefix}.${anc}"
  out="${out_prefix_anc}.pp.txt.gz"
  if [[ ! -f ${out} ]]; then
    Rscript ${rscript} \
      --input_path "${input_path}" \
      --sites "${wes_variants}" \
      --samples "${samples}" \
      --out_prefix "${out_prefix_anc}" \
      --summary_type "long"
  else
    >&2 echo "Skipping ${out}"
  fi
done


