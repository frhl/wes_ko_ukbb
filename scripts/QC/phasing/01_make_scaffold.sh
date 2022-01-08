#!/usr/bin/env bash
#
#$ -N make_scaffold
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/make_scaffold.log
#$ -e logs/make_scaffold.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qc@@short.hga
#$ -t 21

module load foss/2018b
module load HTSlib/1.9-foss-2018b
module load Boost/1.67.0-foss-2018b
module load Rmath/3.5.1-foss-2018b

source utils/bash_utils.sh

readonly in_dir="data/unphased/post-qc"
readonly out_dir="data/mendel"
readonly fam_dir="/well/lindgren/UKBIOBANK/nbaya/resources"

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_wes_200k_scaffold_chr${chr}.vcf.gz"
readonly fam_file="${fam_dir}/ukb11867_pedigree.fam"

readonly scaffold_script="/well/lindgren/flassen/software/phasing/makeScaffold/bin/makeScaffold"

set -x
"${scaffold_script}"\
  --gen ${in_file}\
  --fam ${fam_file}\
  --reg "chr${chr}"\
  --out ${out_file}
set +x

