#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=4G
#SBATCH --time=3-10:00:00
#SBATCH --account=def-barrett

module load singularity
module load nextflow

nextflow run nf-core/methylseq --input "Q08723/guppy/*/*/*_R{1,2}_001.fastq.gz" -profile singularity -with-singularity /home/janayfox/scratch/guppyWGBS/nf-core-methylseq-1.6.1/singularity-images/nfcore-methylseq-1.6.1.img --cytosine_report --fasta /home/janayfox/scratch/guppyWGBS/genome/guppy_genome.fa -resume
