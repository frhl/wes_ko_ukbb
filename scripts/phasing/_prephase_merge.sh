#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh
source utils/hail_utils.sh

module load BCFtools/1.12-GCC-10.3.0

readonly input_list=${1?Error: Missing arg1 (input_list)}
readonly input_type=${2?Error: Missing arg3 (input_type)}
readonly out_prefix=${3?Error: Missing arg3 (output_file)}
readonly out_type=${4?Error: Missing arg4 (output_type)}

readonly out_vcf="${out_prefix}.vcf"
readonly out_vcf_gz="${out_vcf}.gz"
readonly tmp="${out_prefix}.tmp"

rm -f "${out_vcf_gz}"
rm -f "${out_vcf_gz}.tbi"

# remove duplicates (from debugging the functions)
cat ${input_list} | sort | uniq  > ${tmp}

rm_bad_vcf() {
  local _vcf=${1}
  if [ -f ${_vcf} ]; then
    if [ ! -s ${_vcf} ]  || [ $( get_eof_error ${_vcf} ) -gt 0 ]; then
      echo "Removing bad VCF: '${_vcf}' (EOF error or empty file)"
      rm "${_vcf}" "${_vcf}.tbi"
    fi 
 fi
}


# check if VCF is valid
rm_bad_vcf ${out_vcf_gz}


# combine VCFs fast and make tabix
if [ -f "${tmp}" ]; then
  if [ ! -f "${out_vcf}" ]; then
    bcftools merge -l ${tmp} -o "${out_vcf}"
  fi
else
  raise_error "Merge list '${tmp}' does not exist!"
fi

# bgzip
if [ -f "${out_vcf}" ]; then
  if [ ! -f "${out_vcf_gz}" ]; then
    bgzip "${out_vcf}"
  fi
else
  raise_error "File '${out_vcf}' (.vcf) does not exist!"
fi

# index VCF
if [ -f "${out_vcf_gz}" ]; then
  if [ ! -f "${out_vcf_gz}.tbi" ]; then
    make_tabix "${out_vcf_gz}" "tbi"
  fi
else
  raise_error "File '${out_vcf_gz}' (.vcf.gz) does not exist!"
fi

# Combine reads too (assuming eid has already been appended)
readonly file0=$( cat ${tmp} | head -n1 ) 
readonly targetdir=$( dirname ${file0} )
cat ${targetdir}/*.reads > "${out_prefix}.reads"

# clean up temporary files
rm ${tmp}





# remove containing folder and clean up
#readonly prefix_unphased="${out_prefix/_prephased/}"
#readonly unphased_vcf="${prefix_unphased}.vcf.bgz"
#readonly unphased_index="${prefix_unphased}.vcf.bgz.tbi"
#readonly unphased_tmp="${prefix_unphased}.tmp"

# remove directory content and then directory
#rm ${prefix_unphased}/s*
#rmdir ${preifx_unphased}
#rm -f ${unphased_vcf} ${unphased_index} ${unphased_tmp}





