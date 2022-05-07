#!/usr/bin/env bash
#
#$ -N fit_simulated_null
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/fit_simulated_null.log
#$ -e logs/fit_simulated_null.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly chr="${SGE_TASK_ID}"
readonly plink_dir="data/saige/grm/input"
readonly grm_dir="data/saige/grm/input"
readonly covar_dir="data/phenotypes"
readonly pheno_dir="data/simulation/phenotypes"
readonly out_dir="data/simulation/saige/step1"

readonly grm_mtx="${grm_dir}/211102_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly plink_file="${plink_dir}/chunks/ukb_wes_200k_sparse_autosomes_mrg"
readonly covar_file="${covar_dir}/covars2.csv"
readonly covariates=$( cat ${covar_file} )

readonly bash_script="scripts/simulation/_fit_simulated_null.sh"
readonly rscript="utils/saige/step1_fitNULLGLMM.R"

mkdir -p ${out_dir}

fit_phenotypes() {

  local K=0.1
  local h2_nc=${1}
  local h2_co=${2}
  local h2_ko=${3}
  local pi_nc=${4}
  local pi_co=${5}
  local pi_ko=${6}
  local alpha=${7}
  local beta=${8}
  local theta=${9}

  local h2s="${h2_nc}_${h2_co}_${h2_ko}"
  local pis="${pi_nc}_${pi_co}_${pi_ko}"
  local effects="a${alpha}_b${beta}_t${theta}"

  local prefix="ukb_eur_h2_${h2s}_pi_${pis}_K${K}_${effects}_chr${chr}"
  local pheno_file="${pheno_dir}/${prefix}_phenos.tsv.gz"

  local trait_type="quantitative"
  local inv_normalize="TRUE"
  local phenotype="y_cts"
  fit_null
}

fit_null() {
   local out_prefix="${out_dir}/${prefix}_${phenotype}"
   set -x
   qsub -N "_sim_fit_null${SGE_TASK_ID}" \
     -q "short.qe" \
     -pe shmem 1 \
     -t ${tasks} \
     ${bash_script} \
     ${rscript} \
     ${plink_file} \
     ${pheno_file} \
     ${phenotype} \
     ${covariates} \
     ${trait_type} \
     ${inv_normalize} \
     ${out_prefix} \
     ${grm_mtx} \
     ${grm_sam} \
     ${chr}
   set +x
 }

readonly tasks="1-25"

# Absence of CH effects
fit_phenotypes 0.00 0.00 0.00 0.00 0.00 0.00 NA NA NA
fit_phenotypes 0.10 0.00 0.00 0.00 0.00 0.00 NA NA NA
fit_phenotypes 0.00 0.10 0.00 0.00 0.00 0.00 NA NA NA
fit_phenotypes 0.10 0.10 0.00 0.00 0.00 0.00 NA NA NA

# simulate CH effects
#fit_phenotypes 0.00 0.00 0.02 0.00 0.00 0.10 NA NA NA
fit_phenotypes 0.00 0.10 0.02 0.00 0.00 0.10 NA NA NA
fit_phenotypes 0.00 0.10 0.02 0.00 0.00 0.10 NA NA NA

# simulate effects with thetas
fit_phenotypes 0.00 0.00 0.00 0.00 0.00 0.10 NA NA 0.01
fit_phenotypes 0.00 0.00 0.00 0.00 0.00 0.10 NA NA 0.10




