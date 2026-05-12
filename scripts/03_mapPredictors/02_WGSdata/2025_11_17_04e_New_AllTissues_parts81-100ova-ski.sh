#!/bin/bash -l
#SBATCH --job-name=genome_Mapping
#SBATCH --cpus-per-task=56
#SBATCH --account=ypaul1
#SBATCH --output=./tmp/output/genomeRF_o-%A_%a.txt
#SBATCH --error=./tmp/output/genomeRF_e-%A_%a.txt
#SBATCH --array=4,5
#SBATCH --mem=150
#SBATCH --nodelist=beyer-cn08



# Initialize conda
source ~/.bashrc  # Or your conda init file
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R

CHR_IDX=$SLURM_ARRAY_TASK_ID


# Step 2: Process parts 81-100
echo "Processing parts 81-100 for chromosome $CHR_IDX"
Rscript --vanilla '/cellfile/cellnet/MutationModel/scripts/03_mapPredictors/02_WGSdata/2025_11_17_04e_New_AllTissues_parts81-100ova-ski.R' $CHR_IDX "81-100"

echo "Completed processing for chromosome $CHR_IDX"
