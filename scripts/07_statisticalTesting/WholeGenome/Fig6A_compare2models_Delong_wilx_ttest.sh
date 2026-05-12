#!/bin/bash -l
#SBATCH --job-name=Fig6A_stat_tests
#SBATCH --cpus-per-task=8
#SBATCH --account=ypaul1
#SBATCH --output=./tmp/output/Fig1_stat_tests_o-%j.txt
#SBATCH --error=./tmp/output/Fig1_stat_tests_e-%j.txt
#SBATCH --mem=100G
#SBATCH --time=04:00:00

# Initialize conda
source ~/.bashrc
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R


mkdir -p ./tmp/output

Rscript --vanilla 'scripts/07_statisticalTesting/WholeGenome/Fig6A_compare2models_Delong_wilx_ttest.R'