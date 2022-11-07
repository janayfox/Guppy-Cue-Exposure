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

setwd("/scratch/janayfox/guppyWGBS/")

## Prepare tabix files
#create list of file locations
file.list.dev = list("../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC2F4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC2F5.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC2F6.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC3F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC3F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC3M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC3M2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC4F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC4F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC4F3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC4F4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC4F5.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC4M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC5F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC5F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC5F4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC5F5.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC5M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC5M2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC5M3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC5M4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC6F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC6F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC6F3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC6M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC6M2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC6M3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC7F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC7F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC7F3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC7F4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC7F5.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC7F6.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC7M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC7M2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DAC7M3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC2F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC2F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC2M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC2M2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC3F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC3F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC3F3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC3M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC3M2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC3M4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC4F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC4F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC4F3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC4F4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC4F5.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC4M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC4M2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC5F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC5F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC5F3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC5F4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC5F5.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC5M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC5M2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC5M3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC6F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC6F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC6F3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC6F4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC6M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC6M2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC6M4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC7F1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC7F2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC7F3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC7M1.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC7M2.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC7M3.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC7M4.CpG_report.txt.gz",
                 "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/DC7M5.CpG_report.txt.gz")
#create tabix file
myobj=methRead(file.list.dev,
                   sample.id=list("DAC2F4","DAC2F5","DAC2F6","DAC3F1","DAC3F2","DAC3M1",
                                  "DAC3M2","DAC4F1","DAC4F2","DAC4F3","DAC4F4","DAC4F5","DAC4M1",
                                  "DAC5F1", "DAC5F2","DAC5F4","DAC5F5","DAC5M1","DAC5M2","DAC5M3",
                                  "DAC5M4","DAC6F1","DAC6F2","DAC6F3","DAC6M1","DAC6M2",  "DAC6M3",
                                  "DAC7F1","DAC7F2","DAC7F3","DAC7F4","DAC7F5","DAC7F6","DAC7M1",
                                  "DAC7M2","DAC7M3",
                                  "DC2F1","DC2F2", "DC2M1","DC2M2","DC3F1","DC3F2","DC3F3",
                                  "DC3M1","DC3M2","DC3M4","DC4F1","DC4F2","DC4F3","DC4F4",
                                  "DC4F5","DC4M1","DC4M2","DC5F1","DC5F2","DC5F3","DC5F4",
                                  "DC5F5","DC5M1","DC5M2","DC5M3","DC6F1","DC6F2","DC6F3",
                                  "DC6F4","DC6M1","DC6M2","DC6M4","DC7F1","DC7F2","DC7F3",
                                  "DC7M1","DC7M2","DC7M3","DC7M4","DC7M5"),
                   assembly="guppyWGBS_dev_new",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "guppy_dev_DB",
                   mincov = 1
)

#filter out sites in the 99.9th percentile of coverage (PCR bias) and 10x coverage for DMS
myobj.10X=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "10X")

#filter out sites in the 99.9th percentile of coverage (PCR bias) and 3x coverage for DMRs
myobj.3X=filterByCoverage(myobj,lo.count=3,lo.perc=NULL,
                               hi.count=NULL, hi.perc=99.9, suffix = "3X")

#normalize by median coverage
norm.myobj.10X=normalizeCoverage(myobj.10X, method="median")
norm.myobj.3X=normalizeCoverage(myobj.3X, method="median")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = ".backupRData/04_methReadFilter-backup.RData")

q(save="yes")

