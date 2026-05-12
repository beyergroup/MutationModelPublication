#!/bin/bash -l
#SBATCH --job-name=random_forest
#SBATCH --cpus-per-task=56
#SBATCH --account=ypaul1
#SBATCH --output=./tmp/output/RF_o-%A_%a.txt
#SBATCH --error=./tmp/output/RF_e-%A_%a.txt
#SBATCH --array=1-100
#SBATCH --mem=50
#SBATCH --time=20:00:00




# Initialize conda
source ~/.bashrc  # Or your conda init file
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R

# CHR='1'

# Get part number from array task ID
PART=$SLURM_ARRAY_TASK_ID



# Run R script
Rscript --vanilla '/cellfile/cellnet/MutationModel/scripts/04_trainmodels/02_WGSdata/07_wholeGenomeRFPredictionNew.R' $CHR $PART
