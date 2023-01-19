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

setwd("/scratch/janayfox/guppyWGBS_shortterm/24hr/")

## Prepare tabix files
#create lists of file locations
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

#create tabix file
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

#filter out sites in the 99.9th percentile of coverage (PCR bias) 
myobj_24h.filt=filterByCoverage(myobj_24h,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_24h.fem.filt=filterByCoverage(myobj_24h.fem,lo.count=NULL,lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_24h.mal.filt=filterByCoverage(myobj_24h.mal,lo.count=NULL,lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9, suffix = "filt")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/04_shortterm_methReadFilter-backup.RData")

q(save="yes")




