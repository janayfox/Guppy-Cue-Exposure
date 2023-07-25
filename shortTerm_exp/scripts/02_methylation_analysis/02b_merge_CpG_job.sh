#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-01:00
#SBATCH --account=def-barrett

#######################################################
### Goal: Merge coverage on Bismark CpG report files
### Author: Janay Fox
#######################################################

# $1 = Bismark CpG report file 

module load StdEnv/2020
module load python/3.11.2

for i in *.txt.gz; do
    python 02a_merge_CpG.py --input "$i"
done