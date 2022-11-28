#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=agreement
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/agreement.log
#SBATCH --error=logs/agreement.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
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
readonly input_path="data/phased/wes_union_calls/200k/calibration/ukb_shapeit5_whatshap_variants_chr21.PS.txt.gz"

readonly out_dir="data/prephased/wes_union_calls"
readonly out_prefix="${out_dir}/221128_phasing_agreement_n20000"

readonly wes_variants="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

readonly n_samples=20000

mkdir -p ${out_dir}

qsub_agreement() {
  local pp_cutoff=${1}
  local output_path="${out_prefix}_pp${pp_cutoff#*.}.txt"
  qsub -N "_agreement" \
    -o "logs/_agreement.log" \
    -e "logs/_agreement..errors.log" \
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
    ${rscript}
}

qsub_agreement 0.50
qsub_agreement 0.60
qsub_agreement 0.70
qsub_agreement 0.75
qsub_agreement 0.80
qsub_agreement 0.85
qsub_agreement 0.90
qsub_agreement 0.95
qsub_agreement 0.99


