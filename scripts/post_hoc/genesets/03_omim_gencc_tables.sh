#!/usr/bin/env bash
#
#$ -N omim_gencc_tables
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/omim_gencc_tables.log
#$ -e logs/omim_gencc_tables.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/03_omim_gencc_tables.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/230117_evidence"

readonly path_omim="/well/lindgren/flassen/ressources/genesets/genesets/data/omim/morbidmap.txt"
readonly path_gencc="/well/lindgren/flassen/ressources/genesets/genesets/data/genCC/230116_gencc_submissions.tsv"
readonly path_bridge="/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/230117_hgncid_ensembl.txt.gz"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --path_omim ${path_omim} \
  --path_gencc ${path_gencc} \
  --path_bridge ${path_bridge} \
  --out_prefix ${out_prefix}


