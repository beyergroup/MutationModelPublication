#!/bin/bash -l
#SBATCH --cpus-per-task=28
#SBATCH --account=ypaul1
#SBATCH --output=output/exomeRF_o-%j
#SBATCH --error=output/exomeRF_e-%j
#SBATCH --array=1-10

#module unload R-3.5.1
#module load R-4.1.2

# Initialize conda
source ~/.bashrc  # Or your conda init file
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R

Rscript --vanilla --verbose 'scripts/04_trainmodels/01_WESdata/07a_wholeExomeCombined_RFPrediction.R' ${SLURM_ARRAY_TASK_ID}
###this script is run on the new cluster beyer-ln-a
