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

readonly info_dir="data/mt/prefilter/pp90"
readonly info_path="${info_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.frqx.maf.gz"

readonly out_dir="data/mt/prefilter/pp90/encoded"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.recessive.info"

mkdir -p ${out_dir}

export PATH="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/call_chets:${PATH}"

module load BCFtools
#for anno in "pLoF" "damaging_missense" "pLoF_damaging_missense" "synoymous" "other_missense" "non_coding"; do
for anno in pLoF pLoF_damaging_missense; do

  # setup paths 
  out_prefix_anno="${out_prefix}.${anno}"
  map="${out_prefix_anno}.map.txt.gz"
  info="${out_prefix_anno}.info.txt.gz"
  geno="${out_prefix_anno}.geno.txt.gz"
  out="${out_prefix_anno}.encoded.txt.gz"

  # get actual consequence
  if [[ "${anno}" == "pLoF_damaging_missense" ]];then
    zcat ${vep_path} | grep -Ew '(pLoF)|(damaging_missense)' | cut -f3,5 | gzip > ${map}
  else
    zcat ${vep_path} | grep -Ew "${anno}" | cut -f3,5 | gzip > ${map}
  fi

  # get MAC,AN,MAF - note that we seperate by ":" as "|" indicates next haplotype
  zcat ${info_path} | awk -F'\t' '{print $1"\t"$12"\tMAC="$4":AN="$3+$2":AF="$6}' | gzip > ${info}

  # extract genotypes alternate sequences (export by PP)
  zcat ${geno_path} | awk '$4>=0.9 || $4=="" || $4=="."' | gzip > ${geno}

  # call co-occuring variants 

  # additive encoding
  #call_chets --geno ${geno} --gene-map ${map} --gene-collapse "additive" --show-variants | gzip > ${out}
  call_chets --geno ${geno} --gene-map ${map} --info-map ${info} --show-variants | gzip > ${out}

  rm ${map} ${geno} ${info}

done


