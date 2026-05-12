#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --account=cschmal1
#SBATCH --output=output/xgboost_o-%j
#SBATCH --error=output/xgboost_e-%j
#SBATCH --array=8

export HOME=/cellfile/datapublic/cschmalo/cschmal1_home 
conda activate /cellfile/cellnet/MutationModel/conda_envs/r_env
Rscript --vanilla --verbose 'scripts/04_trainmodels/01_WESdata/03b_XGboost.R' ${SLURM_ARRAY_TASK_ID}
