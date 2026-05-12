#!/bin/bash -l
#SBATCH --job-name=Fig5B_stat_tests
#SBATCH --cpus-per-task=14
#SBATCH --account=ypaul1
#SBATCH --output=./tmp/output/Fig5B_stat_tests_o-%j.txt
#SBATCH --error=./tmp/output/Fig5B_stat_tests_e-%j.txt
#SBATCH --mem=100G
#SBATCH --time=04:00:00

mkdir -p ./tmp/output

# Initialize conda
source ~/.bashrc
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R

Rscript --vanilla 'scripts/07_statisticalTesting/WholeExome/Fig5b_compareTumor_universalvsTissueSpecific_Delong_wilx_ttest.R'