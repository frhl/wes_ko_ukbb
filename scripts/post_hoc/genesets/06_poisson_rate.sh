#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=poisson_rate
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/poisson_rate.log
#SBATCH --error=logs/poisson_rate.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --requeue

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/genesets/06_poisson_rate.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/poisson_rate_essential_genesets_review"

readonly in_dir="data/knockouts/tables"
readonly in_file="${in_dir}/combined_annotations_by_sample.new.counts.txt.gz"

readonly geneset_dir="/well/lindgren/flassen/ressources/genesets/genesets/data"
readonly mutation_rates="${geneset_dir}/mutation_rates/samocha2014.txt.gz"
readonly essential_dir="${geneset_dir}/shapeit"
readonly pli="${geneset_dir}/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv"

# get additional (non-mandatory) genesets
readonly omim="${geneset_dir}/omim/230329_morbidmap_by_gene_with_inheritance.txt"
readonly gtex="${geneset_dir}/gtex/GTEx.tstat.tsv"
readonly cancer1="${geneset_dir}/gsea/cancer/gavish_3ca_genes.txt"
readonly cancer2="${geneset_dir}/gsea/cancer/gsea_cgn_genes.txt"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
   --out_prefix "${out_prefix}" \
   --in_count "${in_file}" \
   --file_mutation_rates "${mutation_rates}" \
   --dir_genesets "${essential_dir}" \
   --file_pli ${pli} \
   --file_gtex ${gtex} \
   --file_cancer ${cancer1}

