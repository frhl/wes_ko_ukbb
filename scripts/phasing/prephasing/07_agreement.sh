#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=agreement
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/agreement.log
#SBATCH --error=logs/agreement.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#
#
#$ -N agreement
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/agreement.log
#$ -e logs/agreement.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -V

source utils/bash_utils.sh

readonly rscript="scripts/phasing/prephasing/07_agreement.R"
readonly bashscript="scripts/phasing/prephasing/_agreement.sh"
readonly input_dir="data/prephased/wes_union_calls/phase_conf"
readonly input_path="${input_dir}/ukb_shapeit5_whatshap_chr20.PP.PS.txt.gz" # note: chr21 does not have MAC

readonly out_dir="data/prephased/wes_union_calls/test"
readonly out_prefix="${out_dir}/161222_phasing_agreement_n20000"

readonly wes_variants="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

readonly n_samples=20000

mkdir -p ${out_dir}

qsub_agreement() {
  local pp_cutoff=${1}
  local mac_bin=${2}
  local output_path="${out_prefix}_pp${pp_cutoff#*.}_mac_bin"
  qsub -N "_agreement_${mac_bin}" \
    -o "logs/_agreement.log" \
    -e "logs/_agreement.errors.log" \
    -q "long.qc" \
    -P "lindgren.prjc" \
    -pe shmem 1 \
    -wd $(pwd) \
    ${bashscript} \
    ${input_path} \
    ${output_path} \
    ${n_samples} \
    ${wes_variants} \
    ${pp_cutoff} \
    ${mac_bin} \
    ${rscript}
}

qsub_by_mac_bin() {
  local bin=${1}
  echo "Submitting by mac_bin ${bin}.."
  qsub_agreement 0.50 ${bin}
  qsub_agreement 0.60 ${bin}
  qsub_agreement 0.70 ${bin}
  qsub_agreement 0.75 ${bin}
  qsub_agreement 0.80 ${bin}
  qsub_agreement 0.85 ${bin}
  qsub_agreement 0.90 ${bin}
  qsub_agreement 0.95 ${bin}
  qsub_agreement 0.99 ${bin}
}

qsub_by_mac_bin "singleton"
qsub_by_mac_bin "2-5"
qsub_by_mac_bin "6-10"
qsub_by_mac_bin "11-20"
qsub_by_mac_bin "21-50"
qsub_by_mac_bin "51-100"
qsub_by_mac_bin "101-200"
qsub_by_mac_bin "201-500"
qsub_by_mac_bin "501-1000"
qsub_by_mac_bin "1001-2000"
qsub_by_mac_bin "2001-5000"
qsub_by_mac_bin "5001-10000"



