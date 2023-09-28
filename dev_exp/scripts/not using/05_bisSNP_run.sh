#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=40G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################
### Goal: Run bisSNP
### Author: Janay Fox
#######################################

#need to edit in cpu_cores_number I think? 
module load StdEnv/2020
module load java/17.0.2

#run bisSNP
for i in *.recal.bam; do
    SAMP="${i%%_*}"
    java -Xmx10g -jar BisSNP-1.0.0.jar -R ../genome/guppy_genome.fa -T BisulfiteGenotyper \
    -I "$i" -vfn1 "${SAMP}.cpg.raw.vcf" -vfn2 "${SAMP}.snp.raw.vcf" \
    -stand_call_conf 20 -stand_emit_conf 0 -mmq 30 -mbq 0 -nt 4
done

##double check SAMP thing works ##
#filter fake SNPs
for i in *.snp.raw.vcf; do
    SAMP="${i%%.*}"
    java -Xmx10g -jar BisSNP-1.0.0.jar -R ../genome/guppy_genome.fa -T VCFpostprocess -oldVcf "$i" \
    -newVcf "${i/%.raw.vcf/.filtered.vcf}" -snpVcf "$i" -o "${SAMP}.snp.raw.filter.summary.txt"
done

for i in *.cpg.raw.vcf; do
    SAMP="${i%%.*}"
    java -Xmx10g -jar BisSNP-1.0.0.jar -R ../genome/guppy_genome.fa -T VCFpostprocess -oldVcf "$i" \
    -newVcf "${i/%.raw.vcf/.filtered.vcf}" -snpVcf "${i/%.cpg.raw.vcf/.snp.raw.vcf}" \
    -o "${SAMP}.cpg.raw.filter.summary.txt"
done