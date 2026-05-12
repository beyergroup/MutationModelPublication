#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=cschmal1
#SBATCH --output=output/mapPredictorsDistance_o-%j
#SBATCH --error=output/mapPredictorsDistance_e-%j
#SBATCH --array=1-10

export HOME=/cellfile/datapublic/cschmalo/cschmal1_home 
conda activate /cellfile/cellnet/MutationModel/conda_envs/r_env
Rscript --vanilla --verbose 'scripts/03_mapPredictors/01_WESdata/01b_mapTrainingData_distancePredictors.R' ${SLURM_ARRAY_TASK_ID}
