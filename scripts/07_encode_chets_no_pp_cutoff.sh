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

readonly samples_vcf_dir="data/mt/prefilter/no_pp_cutoff/old"
readonly samples_vcf="${samples_vcf_dir}/ukb_wes_union_calls_200k_chr21.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.vcf.bgz"

readonly in_dir="data/mt/prefilter/no_pp_cutoff/old"
readonly geno_path="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.alt_alleles.txt.gz"

readonly out_dir="data/mt/prefilter/no_pp_cutoff/encoded"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005"

mkdir -p ${out_dir}

export PATH="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/call_chets:${PATH}"

module load BCFtools

readonly samples="${out_prefix}.samples"
bcftools query -l ${samples_vcf} > ${samples}

#for anno in "pLoF_damaging_missense"; do
for anno in "pLoF"; do

  # setup paths 
  out_prefix_anno="${out_prefix}.${anno}"
  map="${out_prefix_anno}.map.txt.gz"
  geno="${out_prefix_anno}.geno.txt.gz"
  out="${out_prefix_anno}.txt.gz"
  vcf="${out_prefix_anno}.vcf.gz"

  # get actual consequence
  if [[ "${anno}" == "pLoF_damaging_missense" ]];then
    zcat ${vep_path} | grep -Ew '(pLoF)|(damaging_missense)' | cut -f3,5 | gzip > ${map}
  else
    zcat ${vep_path} | grep -Ew "${anno}" | cut -f3,5 | gzip > ${map}
  fi

  # extract genotypes alternate sequences (export by PP)
  zcat ${geno_path} | awk '$4>=0.5 || $4=="" || $4=="."' | gzip > ${geno}

  # call variants by haplotypes
  if [[ ! -f ${out} ]]; then
    call_chets --geno ${geno} --gene-map ${map}  --show-variants | gzip > ${out}
  fi
 
  # perform additive encoding
  if [[ ! -f ${vcf} ]]; then
    encode_vcf --input ${out} --min-ac 1 --samples ${samples} --mode "additive" | bgzip > ${vcf}
  fi 

  make_tabix ${vcf} "csi"

  rm ${map} ${geno}

done

rm ${samples}

