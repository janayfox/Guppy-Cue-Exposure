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

setwd("/scratch/janayfox/guppyWGBS_shortterm/1hr/")

## Prepare tabix files
#create lists of file locations
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

#create tabix file
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

#filter out sites in the 99.9th percentile of coverage (PCR bias) 
myobj_1h.filt=filterByCoverage(myobj_1h,lo.count=NULL,lo.perc=NULL,
                               hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_1h.fem.filt=filterByCoverage(myobj_1h.fem,lo.count=NULL,lo.perc=NULL,
                                   hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj_1h.mal.filt=filterByCoverage(myobj_1h.mal,lo.count=NULL,lo.perc=NULL,
                                   hi.count=NULL, hi.perc=99.9, suffix = "filt")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/04_shortterm_methRead_1hr-backup.RData")

q(save="yes")




