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

setwd("/scratch/janayfox/guppyWGBS_shortterm/05hr/")

## Prepare tabix files
#create lists of file locations
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

#create tabix file
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

#filter out sites in the 99.9th percentile of coverage (PCR bias) 
myobj_05h.filt=filterByCoverage(myobj_05h,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_05h.fem.filt=filterByCoverage(myobj_05h.fem,lo.count=NULL,lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_05h.mal.filt=filterByCoverage(myobj_05h.mal,lo.count=NULL,lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9, suffix = "filt")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/04_shortterm_methRead_05hr-backup.RData")

q(save="yes")




