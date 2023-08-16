#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-03:00
#SBATCH --account=def-barrett

#######################################################
### Goal: Add read gruop to BAM files
### Author: Janay Fox
####################################################### 

module load StdEnv/2020
module load picard/2.26.3
module load java/17.0.2
module load samtools/1.17

for i in *L001*; do
    SAMP="${i%%_*}"
    java -Xmx4g -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I="$i" O="${i/%.bam/.withRG.bam}" ID=H3YMTDSX3.1 \
    LB="LIB_${SAMP}" PL=ILLUMINA PU="H3YMTDSX3.1.${SAMP}" SM="$SAMP" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate
done

for i in *L002*; do
    SAMP="${i%%_*}"
    java -Xmx4g -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I="$i" O="${i/%.bam/.withRG.bam}" ID=H3YMTDSX3.2 \
    LB="LIB_${SAMP}" PL=ILLUMINA PU="H3YMTDSX3.2.${SAMP}" SM="$SAMP" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate
done

for i in *withRG.bam; do
    samtools view -H "$i" | grep '^@RG'
done
