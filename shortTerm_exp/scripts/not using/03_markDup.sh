#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=40G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################
### Goal: Mark duplicates
### Author: Janay Fox
#######################################

module load StdEnv/2020
module load picard/2.26.3
module load java/17.0.2
module load samtools/1.17

for i in *withRG*; do
    SAMP="${i%%_*}"
    java -Xmx10g -jar $EBROOTPICARD/picard.jar MarkDuplicates I="$i" O="${i/%.bam/.mdups.bam}" \
    METRICS_FILE="${i/%.bam/mdups.metric.txt}" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
done
