# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

# Note: an error here without any output can indicate that input file is empty!

set -o errexit
set -o nounset

readonly group="original"
readonly samples_dir="/mnt/project/wes_ko_ukbb/data/samples"

for pop in "eur"; do # "afr" "eas" "sas"; do
 in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt/${pop}/${group}"
 out_dir="/wes_ko_ukbb/data/phased/encode_alt/${pop}/${group}"
 dx mkdir -p ${out_dir}
 for mode in "additive" "recessive"; do
   for anno in "pLoF_damaging_missense"; do
      for pp in "0.90" ; do
        for af in "05"; do
          # popmax
          geno="${in_dir}/UKB.wes.merged.phased.qc.${pop}.af${af}.pp${pp}.${group}.${anno}.txt.gz"
          samples="${samples_dir}/UKB.chr21.samples.${pop}.txt"
          females="${samples_dir}/UKB.chrX.samples.${pop}.txt"
          out_prefix_auto="UKB.wes.merged.phased.qc.${pop}.af${af}.pp${pp}.${group}.${anno}.${mode}.auto"
          out_prefix_chrx="UKB.wes.merged.phased.qc.${pop}.af${af}.pp${pp}.${group}.${anno}.${mode}.chrx"
          out="${out_prefix_auto}.vcf.gz"
          echo "GENO = ${geno}"
          echo "SAMPLES = ${samples}"
          echo "OUT = ${out}"
          if [[ $(dx_is_empty $( wo_mnt_project ${geno}) ) -eq "0" ]]; then
            if [[ $(dx_file_exists "${out_dir}/${out}") -eq "0" ]]; then
              dx run app-swiss-army-knife \
                -iimage_file="/docker/call_chets.tar.gz"\
                -icmd="
                  set -x &&
                  zcat ${geno} | grep -wE 'chr[0-9]+' > chr1_22.txt.gz &&
                  zcat ${geno} | grep -w  'chrX' > chrX.txt.gz &&
                  encode_vcf --input chr1_22.txt.gz --min-ac 1 --samples ${samples} --mode ${mode} | bgzip > ${out_prefix_auto}.vcf.gz &&
                  encode_vcf --input chrX.txt.gz --min-ac 1 --samples ${females} --mode ${mode} | bgzip > ${out_prefix_chrx}.vcf.gz &&
                  tabix -C ${out_prefix_auto}.vcf.gz &&
                  tabix -C ${out_prefix_chrx}.vcf.gz &&
                  rm -f chr1_22.txt.gz &&
                  rm -f chrX.txt.gz
                  "\
                --instance-type mem1_ssd1_v2_x8 \
                --folder=".${out_dir}" \
                --priority normal \
                --name encode_${mode}_${anno} -y
             else
               >&2 echo "${out} already exists. Skipping."
             fi
           else
            >&2 echo "File is empty: '${geno}' Exiting loop."  && exit 1
           fi
         done
      done
    done
  done
done


