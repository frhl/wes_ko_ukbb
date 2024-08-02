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

# the following directories and files will be deleted upon successfull exectution
readonly prefix_to_del="${out_prefix/_prephased/}"
readonly dir_to_del="${prefix_to_del}"
readonly vcf_to_del="${dir_to_del}.vcf.bgz"
readonly tbi_to_del="${vcf_to_del}.tbi"
readonly mrg_to_del="${prefix_to_del}.mergelist"
echo "The following directory will be removed: ${dir_to_del}, ${vcf_to_del}, ${mrg_to_del}"

# remove duplicates (from debugging the functions)
cat ${input_list} | sort | uniq  > ${tmp}
readonly n_samples=$( cat ${tmp} | wc )

rm_bad_vcf() {
  local _vcf=${1}
  if [ -f ${_vcf} ]; then
    if [ ! -s ${_vcf} ]  || [ $( get_eof_error ${_vcf} ) -gt 0 ]; then
      echo "Removing bad VCF: '${_vcf}' (EOF error or empty file)"
      rm -f "${_vcf}" "${_vcf}.tbi"
    fi 
 fi
}

validate_samples_in_vcf() {
  local _vcf="${1}"
  local _expt="${2}"
  local _found=$( bcftools query -l "${_vcf}" | wc -l)
  echo "VCF=${_vcf} and expt=${_expt} and found=${_found}"
  if [ "${_expt}" != "${_found}" ]; then
    raise_error "Error: ${_vcf} had ${_found} of ${_expt} expected samples!"
    touch "${_vcf}.ERROR"
    exit 1
  fi
}



# check if VCF is valid
rm_bad_vcf ${out_vcf_gz}

# combine VCFs fast and make tabix
if [ -f "${tmp}" ]; then
  if [ ! -f "${out_vcf_gz}" ]; then
    bcftools merge -l ${tmp} -Oz -o "${out_vcf_gz}"
  fi
else
  raise_error "Merge list '${tmp}' does not exist!"
fi

# ensure that all samples have been merged
validate_samples_in_vcf ${out_vcf_gz} ${n_samples}


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
gzip "${out_prefix}.reads"

# clean up temporary files
echo "Success (${out_prefix})! Deleting temporary files."
rm ${tmp}
rm ${dir_to_del}/s*
rmdir ${dir_to_del}
rm ${vcf_to_del}
rm ${tbi_to_del}
rm ${mrg_to_del}



