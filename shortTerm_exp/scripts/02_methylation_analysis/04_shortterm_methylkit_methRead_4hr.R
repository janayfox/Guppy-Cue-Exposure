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

setwd("/scratch/janayfox/guppyWGBS_shortterm/4hr/")

## Prepare tabix files
#create lists of file locations
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

#create tabix file
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

#filter out sites in the 99.9th percentile of coverage (PCR bias) 
myobj_4h.filt=filterByCoverage(myobj_4h,lo.count=NULL,lo.perc=NULL,
                               hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_4h.fem.filt=filterByCoverage(myobj_4h.mal,lo.count=NULL,lo.perc=NULL,
                                   hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_4h.mal.filt=filterByCoverage(myobj_4h.mal,lo.count=NULL,lo.perc=NULL,
                                   hi.count=NULL, hi.perc=99.9, suffix = "filt")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/04_shortterm_methReadFilter-backup.RData")

q(save="yes")




