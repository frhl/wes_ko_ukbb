#!/usr/bin/env bash
#
# extract genotypes in regions near genes that are significant in primary analysis
#
#$ -N filter_genotypes
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/filter_genotypes.log
#$ -e logs/filter_genotypes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc

module load BCFtools/1.12-GCC-10.3.0

source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly phenotype=${1?Error: Missing arg1 (in dir)}
readonly annotation=${2?Error: Missing arg1 (annotation)}
readonly gene_table=${3?Error: Missing arg1 (in dir)}
readonly final_sample_list=${4?Error: Missing arg1 (in dir)}
readonly out_prefix=${5?Error: Missing arg1 (in dir)}
readonly padding=${6?Error: Missing arg1 (in dir)}
readonly min_maf=${7?Error: Missing arg1 (in dir)}
readonly min_info=${8?Error: Missing arg1 (in dir)}

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/02_filter_genotypes.py"

filter_genotypes() {
  python3 "${hail_script}" \
     --phenotype ${phenotype} \
     --annotation ${annotation} \
     --min_maf ${min_maf} \
     --min_info ${min_info} \
     --padding ${padding} \
     --gene_table ${gene_table} \
     --final_sample_list ${final_sample_list} \
     --out_prefix ${out_prefix}
  print_update "Hail finished writing."
  make_tabix "${out_prefix}.vcf.bgz" "csi"
}

# run analysis
mkdir -p ${out_dir}
set_up_hail
set_up_pythonpath_legacy
filter_genotypes


