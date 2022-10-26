#!/usr/bin/env bash


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly rscript="scripts/phasing/_prephase_sample.R"
readonly hail_script="scripts/phasing/_to_mt.py"

readonly input_path=${1?Error: Missing arg1 (input_path)} 
readonly interval_path=${2?Error: Missing arg1 (input_path)} 
readonly chunk_idx=${3?Error: Missing arg2 (input_path)} 
readonly read_placeholder=${4?Error: Missing arg3 (read_placeholder)} 
readonly mergelist=${5?Error: Missing arg4 (out_prefix)} 
readonly out_prefix=${6?Error: Missing arg4 (out_prefix)} 
readonly sample_idx=$( get_array_task_id ) # one-based index for which sample to phase

readonly refdir="/well/lindgren/flassen/ressources/genome_reference/broad"
readonly grch38="${refdir}/Homo_sapiens_assembly38.fasta"
export REF_PATH="${refdir}/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="${refdir}/ref/cache/%2s/%2s/%s"

# get current sample and read path
set_up_rpy
readonly eid=$( Rscript ${rscript} --chunk_idx ${chunk_idx} --sample_idx ${sample_idx} --interval_path ${interval_path} )
readonly read_path=${read_placeholder/SAMPLE/${eid}}
readonly path_unphased="${out_prefix}/s${sample_idx}_eid${eid}.vcf"
readonly path_unphased_gz="${path_unphased}.gz"

readonly prefix_phased="${out_prefix}/s${sample_idx}_eid${eid}_phased"
readonly path_phased_vcf="${prefix_phased}.vcf"
readonly path_phased_gz="${prefix_phased}.vcf.gz"
readonly path_phased_bgz="${prefix_phased}.vcf.bgz"
readonly out_reads="${out_prefix}/s${sample_idx}_eid${eid}.reads"
readonly tmp_reads="${out_prefix}/s${sample_idx}_eid${eid}.reads.tmp"

# function for extracting a single sample from VCF
extract_sample() {
    local in_vcf=${1}
    local sample=${2}
    local out_vcf=${3}
    if [ ! -f "${out_vcf}.gz" ]; then
      module load BCFtools/1.12-GCC-10.3.0
      bcftools view -s ${sample} -o ${out_vcf} ${in_vcf} && bgzip ${out_vcf}
      make_tabix "${out_vcf}.gz" "tbi"
    fi
}

# extract single sample from VCF
module purge
extract_sample ${input_path} ${eid} ${path_unphased}

if [ ! -f "${path_phased_gz}" ]; then
  SECONDS=0
  module purge
  module load htslib/1.8-gcc5.4.0
  set_up_whatshap
  whatshap phase \
    --reference="${grch38}" \
    --output="${path_phased_vcf}" \
    --output-read-list="${out_reads}" \
    --ignore-read-groups \
    --indels \
    ${path_unphased_gz} \
    ${read_path} \
    && print_update "Finished prephasing ${path_unphased}" ${SECONDS} \
    || raise_error "Error prephasing ${path_unphased}"
  # bgzip sample
  module load BCFtools/1.12-GCC-10.3.0
  bgzip ${path_phased_vcf}
  make_tabix "${path_phased_gz}" "tbi"
fi

# Keep track of files that have been phased
echo "${path_phased_gz}" >> ${mergelist}

