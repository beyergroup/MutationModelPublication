#!/bin/bash -l
#SBATCH --job-name=genome_Mapping
#SBATCH --cpus-per-task=56
#SBATCH --account=ypaul1
#SBATCH --output=./tmp/output/genomeRF_o-%A_%a.txt
#SBATCH --error=./tmp/output/genomeRF_e-%A_%a.txt
#SBATCH --array=3,4,5,6,7,8,9,10,11
#SBATCH --mem=150
#SBATCH --exclude=beyer-cn01,beyer-cn02



# --nodelist=beyer-cn06



# Initialize conda
source ~/.bashrc  # Or your conda init file
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R

CHR_IDX=$SLURM_ARRAY_TASK_ID


# Step 2: Process parts 61-80
echo "Processing parts 61-80 for chromosome $CHR_IDX"
Rscript --vanilla '/cellfile/cellnet/MutationModel/scripts/03_mapPredictors/02_WGSdata/2025_11_17_04d_New_AllTissues_parts61-80.R' $CHR_IDX "61-80"

echo "Completed processing for chromosome $CHR_IDX"
