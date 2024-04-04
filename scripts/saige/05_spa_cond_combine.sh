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

readonly rscript="scripts/saige/05_spa_cond_combine.R"

# phenotypes we are including
readonly header_dir="data/phenotypes"
readonly header_file="${header_dir}/dec22_phenotypes_binary_200k_header.tsv"

# for some reason there is a "bug" in saige, not giving us
# some of the columns when conditioniong on markers, 
# thus we need to map them manually
readonly ref_dir="data/post_hoc/results/2024"
#readonly ref_dir="data/post_hoc/results"
#readonly ref_file="${ref_dir}/176k_sig_saige_sig_prs_pref_N5.txt.gz"
readonly ref_file="${ref_dir}/176k_sig_saige_sig_prs_excl_N5.txt.gz" # ok

# get merged hits that unprocessed
readonly merged_dir="data/conditional/combined/saige"
readonly merged_hits="${merged_dir}/176k_merged_hits_post_cond.txt.gz" # nope

# get merged dominance cond hits that are unprocessed
readonly merged_dominance_dir="data/conditional/dominance/saige"
readonly merged_dominance_hits="${merged_dominance_dir}/176k_merged_hits_post_cond_dominance.txt.gz"

# get merged dominance cond hits that are unprocessed
readonly merged_common_dir="data/conditional/common/saige"
readonly merged_common_hits="${merged_common_dir}/176k_merged_hits_post_common_cond.txt.gz"

# calc bonferroni cutoff
readonly genes_tested=952
readonly phenos_tested=311
readonly bonf_p_cutoff="$(python -c "print(0.05/(${genes_tested}*${phenos_tested}))")"
readonly nom_p_cutoff="$(python -c "print(0.05/(${genes_tested}))")"

# paramters for sutff to include
readonly N_ko_case_cutoff="0"
readonly N_ko_cutoff="5"

readonly out_dir="data/post_hoc/results/2024"
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
  --p_cutoff "${nom_p_cutoff}" \

# write all common hits
out_prefix="${out_dir}/176k_sig_saige_common_cond_sig_pref_prs_combined"
set_up_rpy
Rscript "${rscript}" \
  --ref_file "${ref_file}" \
  --merged_hits "${merged_common_hits}" \
  --N_ko_case_cutoff ${N_ko_case_cutoff} \
  --N_ko_cutoff ${N_ko_cutoff} \
  --out_prefix "${out_prefix}"

out_prefix="${out_dir}/176k_sig_saige_dominance_cond_sig_pref_prs_combined"
set_up_rpy
Rscript "${rscript}" \
  --ref_file "${ref_file}" \
  --merged_hits "${merged_dominance_hits}" \
  --N_ko_case_cutoff ${N_ko_case_cutoff} \
  --N_ko_cutoff ${N_ko_cutoff} \
  --out_prefix "${out_prefix}"


