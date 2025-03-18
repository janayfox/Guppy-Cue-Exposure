#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=60G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca
#SBATCH --job-name=DMS_mal

################################################################################
### Goal: Read alignment files into methylKit and filter to create tabix files of filtered cytosine methylation, call DMS
### Author: Janay Fox
################################################################################

module load StdEnv/2020
module load r/4.3.1

Rscript 10_permutation_DMS_mal.R