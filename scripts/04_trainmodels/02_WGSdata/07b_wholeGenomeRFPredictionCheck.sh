#!/bin/bash -l
#SBATCH --job-name=missing_pos
#SBATCH --cpus-per-task=56
#SBATCH --account=ypaul1
#SBATCH --output=./tmp/output/RF_o-%A_%a.txt
#SBATCH --error=./tmp/output/RF_e-%A_%a.txt
#SBATCH --mem=150G    
#SBATCH --time=100:00:00
#SBATCH --array=1-22

# Define tissues to loop over
# TISSUES=("kidney" "ovary" "prostate" "skin")  # <-- add/remove tissues as needed

TISSUES=("liver")

# Initialize conda
source ~/.bashrc
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R

# SLURM_ARRAY_TASK_ID corresponds to chromosome number (1-22)
CHR=${SLURM_ARRAY_TASK_ID}

# Loop over tissues within each array job
for TISSUE in "${TISSUES[@]}"; do
echo "Running chromosome ${CHR}, tissue ${TISSUE}"
Rscript --vanilla \
'scripts/04_trainmodels/02_WGSdata/07b_wholeGenomeRFPredictionCheck.R' \
${CHR} ${TISSUE}
done
