#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --mem=80g
#SBATCH --output=ASV.out
#SBATCH --error=ASV.error
#SBATCH --job-name=ASV_clustering
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

#load R module
module load foss/2018b
module load R/4.0.0-foss-2018b

#start scripts
./ASV_PacBio.R
