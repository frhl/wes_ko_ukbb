# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"


set -eu
readonly threads=24

#readonly in_dir="/mnt/project//Phasing/PhasingWES/step2_phase_rare/ligated"
readonly in_dir="/mnt/project/wes_ko_ukbb/data/phased/qced/ligated"
readonly out_dir="ews_ko_ukbb/data/phased/export_alt/ligated"
#readonly in_dir="/mnt/project/wes_ko_ukbb/data/phased/qced"
#readonly out_dir="/wes_ko_ukbb/data/phased/export_alt"

dx mkdir -p ${out_dir}
for pop in "eur"; do
  for af in "05"; do
    for CHR in {1..22}; do
      input_vcf="${in_dir}/UKB.wes.chr${CHR}.phased.ligated.qc.${pop}.bcf"
      out_prefix="UKB.wes.chr${CHR}.phased.ligated.qc.${pop}.af${af}"
      if [[ $(dx_file_exists $(wo_mnt_project ${input_vcf})) -eq 1 ]]; then
        if [[ $(dx_file_exists "${out_dir}/${out_prefix}.txt.gz") -eq 0 ]]; then
          dx run app-swiss-army-knife \
            -icmd="
                 bcftools view ${input_vcf} --threads ${threads} --max-af 0.${af} -Ou | bcftools query -i'GT=\"alt\"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT %PP %AC %AN %AF\n]' |  awk '{print \$1\" \"\$2\" \"\$3\" \"\$4\" \"\$5\" \"\$5/\$6}' | gzip > ${out_prefix}.txt.gz && echo 'ok!'
          "\
            --instance-type mem1_ssd1_v2_x8 \
            --folder=".${out_dir}" \
            --priority high \
            --name export_alt_chr${CHR} -y
          else
            >&2 echo "Output '${out_prefix}.tar.gz' already exists. Skipping.."
          fi
        else
          >&2 echo "Input VCF: '${input_vcf}' does not exists!"
        fi
    done
  done
done



