#!/usr/bin/env bash
#
#$ -N clinvar_evidence
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/clinvar_evidence.log
#$ -e logs/clinvar_evidence.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/20_clinvar_evidence.R"
readonly out_dir="derived/tables/clinvar"
readonly out_prefix="${out_dir}/clinvar"

readonly path_sig_hits="data/post_hoc/results/165k_saige_cond_sig_subset_prefer_prs.txt.gz"
readonly path_clinvar="data/knockouts/alt/pp90/clinvar_alleles/ukb_eur_wes_200k_clinvar_chrCHR.txt.gz"
readonly path_phenotype="data/phenotypes/dec22_phenotypes_binary_200k.tsv.gz"
readonly path_info="data/mt/prefilter/final_90_loftee/ukb_wes_union_calls_200k_chrCHR.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.csqs.txt.gz"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --path_hits_to_analyze ${path_sig_hits} \
  --path_phenotype ${path_phenotype} \
  --path_clinvar ${path_clinvar} \
  --path_info ${path_info} \
  --out_prefix ${out_prefix}


