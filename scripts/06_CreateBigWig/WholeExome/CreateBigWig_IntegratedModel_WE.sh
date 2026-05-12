#!/bin/bash -l
#SBATCH --cpus-per-task=28
#SBATCH --account=ypaul1
#SBATCH --output=output/exomeRF_o-%j
#SBATCH --error=output/exomeRF_e-%j


conda activate R 

Rscript --vanilla --verbose 'scripts/06_CreateBigWig/WholeExome/CreateBigWig_IntegratedModel_WE.R'
###this script is run on the new cluster beyer-ln-a
