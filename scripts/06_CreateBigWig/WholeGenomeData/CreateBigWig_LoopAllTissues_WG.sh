#!/bin/bash -l
#SBATCH --job-name=bigwig
#SBATCH --array=0-7
#SBATCH --cpus-per-task=56
#SBATCH --account=ypaul1
#SBATCH --output=./tmp/output/bigwig_o-%A_%a.txt
#SBATCH --error=./tmp/output/bigwig_e-%A_%a.txt
#SBATCH --mem=150
#SBATCH --time=150:00:00

# Map array task ID to tissue name
# TISSUES=(brain breast esophagus kidney liver ovary prostate skin)
TISSUES=(liver)
TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}

echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing tissue: $TISSUE"

# Initialize conda
# Initialize conda
source ~/.bashrc  # Or your conda init file
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R

Rscript --vanilla 'scripts/06_CreateBigWig/WholeGenomeData/CreateBigWig_LoopAllTissues_WG.R' "$TISSUE" 