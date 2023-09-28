#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=40G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################
### Goal: Recalibrate base quality
### Author: Janay Fox
#######################################

#need to edit in cpu_cores_number I think? 
module load StdEnv/2020
module load java/17.0.2

#recalibrate based on three covariates
for i in *withRG*; do
    SAMP="${i%%_*}"
    java -Xmx10g -jar ../../BisSNP-1.0.0.jar -R ../../../genome/guppy_genome.fa -I "$i" \
    -T BisulfiteCountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate \
    -cov CycleCovariate -recalFile "${SAMP}.recalFile_before.csv" -nt 4
done

#write recalibrated base quality score into BAM files
for i in *withRG*; do
    SAMP="${i%%_*}"
    java -Xmx10g -jar ../../BisSNP-1.0.0.jar -R ../genome/guppy_genome.fa -I "$i" \
    -o "${i/%.bam/.recal.bam}" -T BisulfiteTableRecalibration \
    -recalFile "${SAMP}.recalFile_before.csv" -maxQ 40
done

#validate recalibration step 
# for i in *.recal.bam; do
#     SAMP="${i%%_*}"
#     java -Xmx10g -jar ../../BisSNP-1.0.0.jar -R ../genome/guppy_genome.fa -I "$i" \
#     -T BisulfiteCountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate \
#     -cov CycleCovariate -recalFile "${SAMP}.recalFile_after.csv" -nt 4
# done