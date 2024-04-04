#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J coxph_models
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/logs/coxph_models-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

# If using bulk submission via R script:
echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}



module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

#source /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/bash_utils.sh
#set_up_rpy


Rscript /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/scripts/survival/run_original/updated_2311_subset_run_coxph_models.R \
--geneNames=${geneNames} \
--chunk=${chunk} \
--fileType=${fileType} \
--refGroup=${refGroup}

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0


