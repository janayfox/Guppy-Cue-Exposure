#####################################################################################################################
### Goal: Read alignment files into methylKit and filter to create tabix files of filtered cytosine methylation
### Author: Janay Fox
### R script
#####################################################################################################################

## Set up  ##
#install packages 
#install.packages("S4Vectors")
#install.packages("IRanges")
#install.packages("GenomicRanges")
#install.packages("methylKit")

#load packages 
library("S4Vectors", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("GenomicRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")
library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS_shortterm/")

## Prepare tabix files
#create lists of file locations
file.list.all = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC10F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC10M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC11F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC11M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC13F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC13M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC14F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC14M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC15F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC15M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC16F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC16M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC1F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC1M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC2F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC2M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC3F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC3M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC4F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC4M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC5F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC5M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC6F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC6M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC7F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC7M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC8F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC8M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC9F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC9M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C10F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C10M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C11F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C11M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C12F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C12M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C13F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C13M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C14F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C14M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C15F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C15M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C1F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C1M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C2F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C2M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C3F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C3M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C4F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C4M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C5F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C5M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C6F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C6M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C7F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C7M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C8F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C8M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C9F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C9M.CpG_report.txt.gz")

file.list.all.fem = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC10F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC11F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC13F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC14F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC15F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC16F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC1F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC2F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC3F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC4F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC5F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC6F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC7F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC8F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC9F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C10F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C11F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C12F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C13F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C14F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C15F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C1F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C2F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C3F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C4F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C5F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C6F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C7F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C8F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C9F.CpG_report.txt.gz")

file.list.all.mal = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC10M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC11M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC13M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC14M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC15M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC16M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC1M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC2M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC3M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC4M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC5M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC6M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC7M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC8M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC9M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C10M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C11M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C12M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C13M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C14M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C15M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C1M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C2M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C3M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C4M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C5M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C6M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C7M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C8M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C9M.CpG_report.txt.gz")

file.list.05h = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC10F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC10M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC16F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC16M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC3F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC3M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C10F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C10M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C13F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C13M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C2F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C2M.CpG_report.txt.gz")

file.list.05h.fem = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC10F.CpG_report.txt.gz",
                      "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC16F.CpG_report.txt.gz",
                      "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC3F.CpG_report.txt.gz",
                      "../results/bismark_methylation_calls/stranded_CpG_report/ST2C10F.CpG_report.txt.gz",
                      "../results/bismark_methylation_calls/stranded_CpG_report/ST2C13F.CpG_report.txt.gz",
                      "../results/bismark_methylation_calls/stranded_CpG_report/ST2C2F.CpG_report.txt.gz")

file.list.05h.mal = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC10M.CpG_report.txt.gz",
                      "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC16M.CpG_report.txt.gz",
                      "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC3M.CpG_report.txt.gz",
                      "../results/bismark_methylation_calls/stranded_CpG_report/ST2C10M.CpG_report.txt.gz",
                      "../results/bismark_methylation_calls/stranded_CpG_report/ST2C13M.CpG_report.txt.gz",
                      "../results/bismark_methylation_calls/stranded_CpG_report/ST2C2M.CpG_report.txt.gz")

file.list.1h = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC14F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC14M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC2F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC2M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC9F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC9M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C14F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C14M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C3F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C3M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C9F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C9M.CpG_report.txt.gz")

file.list.1h.fem = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC14F.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC2F.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC9F.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C14F.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C3F.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C9F.CpG_report.txt.gz")

file.list.1h.mal = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC14M.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC2M.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC9M.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C14M.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C3M.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C9M.CpG_report.txt.gz")

file.list.4h = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC13F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC13M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC1F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC1M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC7F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC7M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C12F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C12M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C1F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C1M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C7F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C7M.CpG_report.txt.gz")

file.list.4h.fem = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC13F.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC1F.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC7F.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C12F.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C1F.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C7F.CpG_report.txt.gz")

file.list.4h.mal = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC13M.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC1M.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC7M.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C12M.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C1M.CpG_report.txt.gz",
                    "../results/bismark_methylation_calls/stranded_CpG_report/ST2C7M.CpG_report.txt.gz")


file.list.24h = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC15F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC15M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC4F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC4M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC8F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC8M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C15F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C15M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C4F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C4M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C8F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C8M.CpG_report.txt.gz")

file.list.24h.fem = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC15F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC4F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC8F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C15F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C4F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C8F.CpG_report.txt.gz")

file.list.24h.mal = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC15M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC4M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC8M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C15M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C4M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C8M.CpG_report.txt.gz")

file.list.72h = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC11F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC11M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC5F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC5M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC6F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC6M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C11F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C11M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C5F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C5M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C6F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C6M.CpG_report.txt.gz")

file.list.72h.fem = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC11F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC5F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC6F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C11F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C5F.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C6F.CpG_report.txt.gz")

file.list.72h.mal = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC11M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC5M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2AC6M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C11M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C5M.CpG_report.txt.gz",
                     "../results/bismark_methylation_calls/stranded_CpG_report/ST2C6M.CpG_report.txt.gz")

#create tabix file
myobj_all=methRead(file.list.all,
                   sample.id=list("ST2AC10F","ST2AC10M","ST2AC11F","ST2AC11M","ST2AC13F",
                                  "ST2AC13M","ST2AC14F","ST2AC14M","ST2AC15F","ST2AC15M",
                                  "ST2AC16F","ST2AC16M","ST2AC1F","ST2AC1M","ST2AC2F",
                                  "ST2AC2M","ST2AC3F","ST2AC3M","ST2AC4F","ST2AC4M",
                                  "ST2AC5F","ST2AC5M","ST2AC6F","ST2AC6M","ST2AC7F",
                                  "ST2AC7M","ST2AC8F","ST2AC8M","ST2AC9F","ST2AC9M",
                                  "ST2C10F","ST2C10M","ST2C11F","ST2C11M","ST2C12F",
                                  "ST2C12M","ST2C13F","ST2C13M","ST2C14F","ST2C14M",
                                  "ST2C15F","ST2C15M","ST2C1F","ST2C1M","ST2C2F",
                                  "ST2C2M","ST2C3F","ST2C3M","ST2C4F","ST2C4M",
                                  "ST2C5F","ST2C5M","ST2C6F","ST2C6M","ST2C7F",
                                  "ST2C7M","ST2C8F","ST2C8M","ST2C9F", "ST2C9M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_all_DB"
)

myobj_all.fem=methRead(file.list.all.fem,
                   sample.id=list("ST2AC10F","ST2AC11F","ST2AC13F","ST2AC14F","ST2AC15F",
                                  "ST2AC16F","ST2AC1F","ST2AC2F","ST2AC3F","ST2AC4F",
                                  "ST2AC5F","ST2AC6F","ST2AC7F","ST2AC8F","ST2AC9F",
                                  "ST2C10F","ST2C11F","ST2C12F","ST2C13F","ST2C14F",
                                  "ST2C15F","ST2C1F","ST2C2F","ST2C3F","ST2C4F",
                                  "ST2C5F","ST2C6F","ST2C7F","ST2C8F","ST2C9F"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_allfem_DB"
)

myobj_all.mal=methRead(file.list.all.mal,
                   sample.id=list("ST2AC10M","ST2AC11M","ST2AC13M","ST2AC14M","ST2AC15M",
                                  "ST2AC16M","ST2AC1M","ST2AC2M","ST2AC3M","ST2AC4M",
                                  "ST2AC5M","ST2AC6M","ST2AC7M","ST2AC8M","ST2AC9M",
                                  "ST2C10M","ST2C11M","ST2C12M","ST2C13M","ST2C14M",
                                  "ST2C15M","ST2C1M","ST2C2M","ST2C3M","ST2C4M",
                                  "ST2C5M","ST2C6M","ST2C7M","ST2C8M", "ST2C9M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_allmal_DB"
)

myobj_05h=methRead(file.list.05h,
                   sample.id=list("ST2AC10F","ST2AC10M","ST2AC16F","ST2AC16M","ST2AC3F","ST2AC3M",
                                  "ST2C10F","ST2C10M","ST2C13F","ST2C13M","ST2C3F","ST2C3M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,1,1,1,
                               0,0,0,0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_05h_DB"
)

myobj_05h.fem=methRead(file.list.05h.fem,
                   sample.id=list("ST2AC10F","ST2AC16F","ST2AC3F",
                                  "ST2C10F","ST2C13F","ST2C3F"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,
                               0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_05hfem_DB"
)

myobj_05h.mal=methRead(file.list.05h.mal,
                   sample.id=list("ST2AC10M","ST2AC16M","ST2AC3M",
                                  "ST2C10M","ST2C13M","ST2C3M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,
                               0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_05hmal_DB"
)

myobj_1h=methRead(file.list.1h,
                   sample.id=list("ST2AC14F","ST2AC14M","ST2AC2F","ST2AC2M","ST2AC9F","ST2AC9M",
                                  "ST2C14F","ST2C14M","ST2C3F","ST2C3M","ST2C9F", "ST2C9M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,1,1,1,
                               0,0,0,0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_1h_DB"
)

myobj_1h.fem=methRead(file.list.1h.fem,
                  sample.id=list("ST2AC14F","ST2AC2F","ST2AC9F",
                                 "ST2C14F","ST2C3F","ST2C9F"),
                  assembly="guppyWGBS_shortterm",
                  pipeline="bismarkCytosineReport",
                  treatment=c(1,1,1,
                              0,0,0),
                  context="CpG",
                  mincov = 3,
                  dbtype = "tabix",
                  dbdir = "shortterm_1hfem_DB"
)

myobj_1h.mal=methRead(file.list.1h.mal,
                  sample.id=list("ST2AC14M","ST2AC2M","ST2AC9M",
                                 "ST2C14M","ST2C3M","ST2C9M"),
                  assembly="guppyWGBS_shortterm",
                  pipeline="bismarkCytosineReport",
                  treatment=c(1,1,1,
                              0,0,0),
                  context="CpG",
                  mincov = 3,
                  dbtype = "tabix",
                  dbdir = "shortterm_1hmal_DB"
)

myobj_4h=methRead(file.list.4h,
                   sample.id=list("ST2AC13F","ST2AC13M","ST2AC1F","ST2AC1M","ST2AC7F","ST2AC7M",
                                  "ST2C12F","ST2C12M","ST2C1F","ST2C1M","ST2C7F","ST2C7M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,1,1,1,
                               0,0,0,0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_4h_DB"
)

myobj_4h.fem=methRead(file.list.4h.fem,
                  sample.id=list("ST2AC13F","ST2AC1F","ST2AC7F",
                                 "ST2C12F","ST2C1F","ST2C7F"),
                  assembly="guppyWGBS_shortterm",
                  pipeline="bismarkCytosineReport",
                  treatment=c(1,1,1,
                              0,0,0),
                  context="CpG",
                  mincov = 3,
                  dbtype = "tabix",
                  dbdir = "shortterm_4hfem_DB"
)

myobj_4h.mal=methRead(file.list.4h.mal,
                  sample.id=list("ST2AC13M","ST2AC1M","ST2AC7M",
                                 "ST2C12M","ST2C1M","ST2C7M"),
                  assembly="guppyWGBS_shortterm",
                  pipeline="bismarkCytosineReport",
                  treatment=c(1,1,1,
                              0,0,0),
                  context="CpG",
                  mincov = 3,
                  dbtype = "tabix",
                  dbdir = "shortterm_4hmal_DB"
)

myobj_24h=methRead(file.list.24h,
                   sample.id=list("ST2AC15F","ST2AC15M","ST2AC4F","ST2AC4M","ST2AC8F","ST2AC8M",
                                  "ST2C15F","ST2C15M","ST2C4F","ST2C4M","ST2C8F","ST2C8M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,1,1,1,
                               0,0,0,0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_24h_DB"
)

myobj_24h.fem=methRead(file.list.24h.fem,
                   sample.id=list("ST2AC15F","ST2AC4F","ST2AC8F",
                                  "ST2C15F","ST2C4F","ST2C8F"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,
                               0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_24hfem_DB"
)

myobj_24h.mal=methRead(file.list.24h.mal,
                   sample.id=list("ST2AC15M","ST2AC4M","ST2AC8M",
                                  "ST2C15M","ST2C4M","ST2C8M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,
                               0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_24hmal_DB"
)

myobj_72h=methRead(file.list.72h,
                   sample.id=list("ST2AC11F","ST2AC11M","ST2AC5F","ST2AC5M","ST2AC6F","ST2AC6M",
                                  "ST2C11F","ST2C11M","ST2C5F","ST2C5M","ST2C6F","ST2C6M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,1,1,1,
                               0,0,0,0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_72h_DB"
)

myobj_72h.fem=methRead(file.list.72h.fem,
                   sample.id=list("ST2AC11F","ST2AC5F","ST2AC6F",
                                  "ST2C11F","ST2C5F","ST2C6F"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,
                               0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_72hfem_DB"
)

myobj_72h.mal=methRead(file.list.72h.mal,
                   sample.id=list("ST2AC11M","ST2AC5M","ST2AC6M",
                                  "ST2C11M","ST2C5M","ST2C6M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,
                               0,0,0),
                   context="CpG",
                   mincov = 3,
                   dbtype = "tabix",
                   dbdir = "shortterm_72hmal_DB"
)


#filter out sites in the 99.9th percentile of coverage (PCR bias) 
myobj_all.filt=filterByCoverage(myobj_all,lo.count=NULL,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_all.fem.filt=filterByCoverage(myobj_all.fem,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_all.mal.filt=filterByCoverage(myobj_all.mal,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_05h.filt=filterByCoverage(myobj_05h,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_05h.fem.filt=filterByCoverage(myobj_05h.fem,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_05h.mal.filt=filterByCoverage(myobj_05h.mal,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_1h.filt=filterByCoverage(myobj_1h,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_1h.fem.filt=filterByCoverage(myobj_1h.fem,lo.count=NULL,lo.perc=NULL,
                               hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_1h.mal.filt=filterByCoverage(myobj_1h.mal,lo.count=NULL,lo.perc=NULL,
                               hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_4h.filt=filterByCoverage(myobj_4h,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_4h.fem.filt=filterByCoverage(myobj_4h.mal,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_4h.mal.filt=filterByCoverage(myobj_4h.mal,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_24h.filt=filterByCoverage(myobj_24h,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_24h.fem.filt=filterByCoverage(myobj_24h.fem,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_24h.mal.filt=filterByCoverage(myobj_24h.mal,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_72h.filt=filterByCoverage(myobj_72h,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_72h.fem.filt=filterByCoverage(myobj_72h.fem,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_72h.mal.filt=filterByCoverage(myobj_72h.mal,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/04_shortterm_methReadFilter-backup.RData")

q(save="yes")




