#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=45G
#SBATCH --time=1-00:00
#SBATCH --account=rrg-barrett
#SBATCH --job-name=methRead
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load StdEnv/2020
module load r/4.3.1

Rscript 04a_shortterm_methylkit_methReadCallDMS_1hr.R