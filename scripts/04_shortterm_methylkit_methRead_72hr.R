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

setwd("/scratch/janayfox/guppyWGBS_shortterm/72hr/")

## Prepare tabix files
#create lists of file locations
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




