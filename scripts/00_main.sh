# This file contains the ordered workflow of the scripts used in this project

DATE="210224"

# 1) run VEP on 200k WES data to generate variant effect annotations.
RUN1="vep"
qsub -N ${RUN1} \
     -o "${DATE}.${RUN1}.log" \
     -e "${DATE}.${RUN1}.errors.log" \
     -q short.qc \
     -pe shmem 8 \
     -t 1-22 01_vep.sh

# 2) extract the VEP variants from bed files in one
RUN2="combine_vep"
qsub -hold_jid ${RUN1} \
     -N ${RUN2} \
     -o "${DATE}.${RUN2}.log" \
     -e "${DATE}.${RUN2}.errors.log" \
     -q short.qc \
     -pe shmem 1 \
     02_combine_vep.sh


