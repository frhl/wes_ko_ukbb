#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=encode_chets
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/encode_chets.log
#SBATCH --error=logs/encode_chets.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22
#SBATCH --dependency="afterok:40456440"

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/vcf_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly vep_version="95"
readonly vep_dir="data/vep/vep${vep_version}/worst_csqs"
readonly vep_path="${vep_dir}/UKB.chr${chr}.exome_array.variants_only.vep${vep_version}.csqs.worst_csq_by_gene_canonical.original.txt.gz"

readonly in_dir="data/mt/prefilter/pp90"
readonly geno_path="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.alt_alleles.txt.gz"

readonly out_dir="data/mt/prefilter/pp90/encoded"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt"

mkdir -p ${out_dir}

export PATH="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/call_chets:${PATH}"

module load BCFtools
for anno in "pLoF" "damaging_missense" "pLoF_damaging_missense" "synoymous" "other_missense" "non_coding"; do

  # setup paths 
  out_prefix_anno="${out_prefix}.${anno}"
  map="${out_prefix_anno}.map.txt.gz"
  geno="${out_prefix_anno}.geno.txt.gz"
  out="${out_prefix_anno}.encoded.txt.gz"

  # get actual consequence
  if [[ "${anno}" == "pLoF_damaging_missense" ]];then
    zcat ${vep_path} | grep -Ew '(pLoF)|(damaging_missense)' | cut -f3,5 | gzip > ${map}
  else
    zcat ${vep_path} | grep -Ew "${anno}" | cut -f3,5 | gzip > ${map}
  fi

  # extract genotypes alternate sequences (export by PP)
  zcat ${geno_path} | awk '$4>=0.5 || $4=="" || $4=="."' | gzip > ${geno}

  # call co-occuring variants 
  call_chets --geno ${geno} --gene-map ${map} --show-variants | gzip > ${out}

  rm ${map} ${geno}

done


