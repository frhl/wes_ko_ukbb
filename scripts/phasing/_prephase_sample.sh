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

set_up_rpy
# get current sample and read path
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
  # append sample ID to read table so that we can follow who is wo

fi

#cat ${out_reads} | awk -v eid="${eid}" 'NR==1{$9="eid";print;next} {$9=eid}1' | column -t > ${tmp_reads}

echo "${path_phased_gz}" >> ${mergelist}

#rm -f "${path_phased_gz}"
#rm -f "${path_phased_gz}.tbi"

#if [ ! -f "${path_phased_bgz}" ]; then
#set +eu
#module purge
#set_up_hail
#set_up_pythonpath_legacy
#python3 ${hail_script} \
#  --input_path "${prefix_phased}.mt" \
#  --input_type "mt" \
#  --out_prefix ${prefix_phased} \
#  --out_type "vcf"
#set -eu
#fi







#module load BCFtools/1.12-GCC-10.3.0
#bgzip -d ${path_phased_gz}

#rm -f "${path_phased_gz}.tbi"
#rm -f "${path_phased_gz}"
#rm -f "${path_phased_bgz}"


#bgzip ${path_phased}
#make_tabix "${path_phased_gz}" "tbi"






# uncompress
#zcat ${path_phased_gz} > ${path_phased}
#bgzip -d ${path_phased_gz}

#gunzip -c ${path_phased} | bgzip > ${path_phased_bgz}


#make_tabix "${path_phased_bgz}" "tbi"

#bgzip ${path_phased}

# make tabix
#if [ ! -f "${path_phased_gz}.tbi" ]; then
#  make_tabix "${path_phased_gz}" "tbi"
#fi

# tidy up
#rm "${path_unphased_gz}" "${path_unphased_gz}.tbi"





