#!/usr/bin/env bash
#
# Get a nice table of the final hits after conditioning on common and rare variants (and PRS)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=spa_cond_combine
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/spa_cond_combine.log
#SBATCH --error=logs/spa_cond_combine.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/13_spa_cond_combine.R"

# phenotypes we are including
readonly header_dir="data/phenotypes"
readonly header_file="${header_dir}/dec22_phenotypes_binary_200k_header.tsv"

# for some reason there is a "bug" in saige, not giving us
# some of the columns when conditioniong on markers, 
# thus we need to map them manually
readonly ref_dir="data/post_hoc/results"
readonly ref_file="${ref_dir}/176k_sig_saige_sig_prs_pref.txt.gz"

# get merged hits that unprocessed
readonly merged_dir="data/conditional/combined/saige"
readonly merged_hits="${merged_dir}/176k_merged_hits_post_cond.txt.gz"

# calc bonferroni cutoff
readonly genes_tested=958
readonly phenos_tested=311
readonly p_cutoff="$(python -c "print(0.05/(${genes_tested}*${phenos_tested}))")"

# paramters for sutff to include
readonly N_ko_case_cutoff="2"
readonly N_ko_cutoff="5"

readonly out_dir="data/post_hoc/results"
mkdir -p ${out_dir}

# write all hits interrogated regardless of final P
out_prefix="${out_dir}/176k_sig_saige_cond_all_pref_prs_combined"
set_up_rpy
Rscript "${rscript}" \
  --ref_file "${ref_file}" \
  --merged_hits "${merged_hits}" \
  --N_ko_case_cutoff ${N_ko_case_cutoff} \
  --N_ko_cutoff ${N_ko_cutoff} \
  --out_prefix "${out_prefix}"

# only write hits with sig final P
out_prefix="${out_dir}/176k_sig_saige_cond_sig_pref_prs_combined"
set_up_rpy
Rscript "${rscript}" \
  --ref_file "${ref_file}" \
  --merged_hits "${merged_hits}" \
  --N_ko_case_cutoff ${N_ko_case_cutoff} \
  --N_ko_cutoff ${N_ko_cutoff} \
  --out_prefix "${out_prefix}" \
  --p_cutoff ${p_cutoff} \

