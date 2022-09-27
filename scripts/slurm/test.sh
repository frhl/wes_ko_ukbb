#!/usr/bin/env bash
#
#BATCH -A lindgren.prj
#SBATCH -J slurm_test
#SBATCH -D /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/slurm_test.log
#SBATCH --error=logs/slurm_test.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=20-21
#SBATCH --requeue

echo -e "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"
echo -e "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo -e "SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}"
echo "SLURM_JOB_PARTITION=${SLURM_JOB_PARTITION}"
echo "SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}"

source utils/qsub_utils.sh
source utils/hail_utils.sh

set_up_hail
set_up_pythonpath_legacy

readonly hail_script="scripts/slurm/test.py"
readonly message="hello"
readonly message2="world"
>&2 echo "the current array is ${SLURM_ARRAY_TASK_ID}"


#python3 ${hail_script} --arg1 ""hello!""   # <-- error (illegal instruction)
#python3 ${hail_script} --arg1 "${message}" # seem to be working
#python3 ${hail_script} --arg1 "my_msg_without_spaces_is_${message}" # seem to be working
#python3 ${hail_script} --arg1 "my msg with spaces ${message}" # seem to be working??
#python3 ${hail_script} --arg1 "my msg with some weird character/${message}" # seem to be working??
#python3 ${hail_script} --arg1 ${message} --arg2 ${message2}
SECONDS=0
python3 "${hail_script}" \
  --arg1 "${message}" \
  --arg2 ${message2}
