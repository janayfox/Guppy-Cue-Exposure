#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-10:00:00
#SBATCH --account=def-barrett

tar -zcvf dev.tar.gz ./SRA/dev 
tar -zcvf st.tar.gz ./SRA/st 