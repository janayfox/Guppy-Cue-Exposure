#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-03:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Merge BAM files and run BS-SNPer
### Author: Janay Fox
####################################################### 


module load StdEnv/2020
module load samtools/1.17
module load perl/5.30.2 
module load gcc/11.3.0

## Sort BAM files
for i in *.bam; do
   samtools stats "$i" | grep "is sorted:"
done

for i in *.bam; do
   samtools sort "$i" > "${i/%.bam/.sorted.bam}"
done

##Merge BAM files 
samtools merge shortterm_merged_sorted_deduplicated.bam *.bam

#view header to check 
samtools view shortterm_merged_sorted_deduplicated.bam | head

#create index file 
samtools index -b shortterm_merged_sorted_deduplicated.bam

##Run BS-SNPer with default settings 
perl /home/janayfox/scratch/guppyWGBS/BS-SNPer/BS-Snper-master/BS-Snper.pl \
/home/janayfox/scratch/guppyWGBS/BS-SNPer/shortterm_merged_sorted_deduplicated.bam \
--fa /home/janayfox/scratch/guppyWGBS/genome/guppy_genome.fa \
--output /home/janayfox/scratch/guppyWGBS/BS-SNPer/shortterm_SNP_candidates.txt \
--methcg /home/janayfox/scratch/guppyWGBS/BS-SNPer/shortterm_CpG_meth_info.tab \
--methchg /home/janayfox/scratch/guppyWGBS/BS-SNPer/shortterm_CHG_meth_info.tab \
--methchh /home/janayfox/scratch/guppyWGBS/BS-SNPer/shortterm_CHH_meth_info.tab \
--minhetferq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 \
--minread2 2 --errorate 0.02 --mapvalue 20 >shortterm_SNP-results.vcf 2>shortterm_ERR.log

#extract only C to T SNPs
grep 'C'$'\t''T' shortterm_SNP-results.vcf >shortterm_CT_SNP.tab