#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-03:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

################################################################################
### Goal: Run gene annotation scripts
### Author: Janay Fox
################################################################################

module load StdEnv/2020
module load r/4.3.1

Rscript 06a_geneAnnotation.R