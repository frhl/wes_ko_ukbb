#!/usr/bin/env bash
#
# extract genotypes in regions near genes that are significant in primary analysis
#
#$ -N _filter_genotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_filter_genotypes.log
#$ -e logs/_filter_genotypes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc

module load BCFtools/1.12-GCC-10.3.0

source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly gene_table=${1?Error: Missing arg1 (Gene table intervals)}
readonly final_sample_list=${2?Error: Missing arg2 (samples to include)}
readonly out_prefix=${3?Error: Missing arg3 (out prefix)}
readonly padding=${4?Error: Missing arg4 (padding to add upstream/downstream of gens)}
readonly min_maf=${5?Error: Missing arg5 (Filter variants by min MAF)}
readonly min_info=${6?Error: Missing arg1 (Filter variants by min INFO)}
readonly pheno_file=${7?Error: Missing arg7 (pheno_file)}
readonly trait=${8?Error: Missing arg8 (trait)}
readonly phenotype=${9?Error: Missing arg9 (phenotype)}

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/common/02_filter_genotypes.py"

filter_genotypes() 
{
  if [ ! -f ${out_prefix}.vcf.bgz ]; then
    python3 "${hail_script}" \
       --min_maf ${min_maf} \
       --min_info ${min_info} \
       --padding ${padding} \
       --gene_table ${gene_table} \
       --extract ${final_sample_list} \
       --out_prefix ${out_prefix} \
       --pheno_file ${pheno_file} \
       --phenotype ${phenotype} \
       --trait ${trait} \
       --min_maf_by_case_control \
       && print_update "Finished filtering imputed genotypes ${out_prefix}" ${SECONDS} \
       || raise_error "Filtering imputed genotypes for for ${out_prefix} failed!"
  else
    >%2 echo "${out_prefix}.vcf.bgz already exists. Skipping.."
  fi
}

# run analysis
set_up_hail
set_up_pythonpath_legacy
filter_genotypes
make_tabix "${out_prefix}.vcf.bgz" "csi"

