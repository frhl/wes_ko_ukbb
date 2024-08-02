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

readonly rscript="scripts/post_hoc/06_clinvar_evidence.R"
readonly path_clinvar="data/knockouts/alt/pp90/clinvar_alleles/ukb_eur_wes_200k_clinvar_chrCHR.txt.gz"
readonly path_info="data/mt/prefilter/final_90_loftee/ukb_wes_union_calls_200k_chrCHR.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.csqs.txt.gz"

readonly path_sig_hits="data/post_hoc/results/176k_saige_cond_sig_subset_prefer_prs.txt.gz"
readonly path_phenotype="data/phenotypes/dec22_phenotypes_binary_200k.tsv.gz"

readonly path_tte_sig_hits="data/survival/results/2302_sig_assoc.txt.gz"
readonly path_tte_phenotype="data/phenotypes/jan23_tte_spiro_phenotypes_logical.txt.gz"

readonly out_dir="derived/tables/clinvar"
readonly out_tte_prefix="${out_dir}/176k_clinvar_tte"
readonly out_prefix="${out_dir}/176k_clinvar"


mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --path_hits_to_analyze ${path_tte_sig_hits} \
  --path_phenotype ${path_tte_phenotype} \
  --path_clinvar ${path_clinvar} \
  --path_info ${path_info} \
  --out_prefix ${out_tte_prefix} \
  --path_hits_type "hazard"


#set_up_rpy
#Rscript "${rscript}" \
#  --path_hits_to_analyze ${path_sig_hits} \
#  --path_phenotype ${path_phenotype} \
#  --path_clinvar ${path_clinvar} \
#  --path_info ${path_info} \
#  --out_prefix ${out_prefix}


