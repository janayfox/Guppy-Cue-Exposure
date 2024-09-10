###########################################################################################
### Goal: Run DMR analysis
### Author: Janay Fox
### R script
###########################################################################################

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
library("data.table", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/methylKit/dev/perc20")

## Prepare tabix files
#create list of file locations
file.list.dev.all = list("../../../mergedCov/dev/DAC2F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC2F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC2F6.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC3F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC3F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC3M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC3M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC4F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC4F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC4F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC4F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC4F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5M3.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5M4.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6M3.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F6.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7M3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC2F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC2F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC2M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC2M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3M4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5M3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6M4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7M3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7M4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7M5.CpG_merged.cov")

file.list.dev.fem = list("../../../mergedCov/dev/DAC2F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC2F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC2F6.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC3F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC3F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC4F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC4F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC4F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC4F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC4F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7F6.CpG_merged.cov",
                     "../../../mergedCov/dev/DC2F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC2F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5F5.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6F3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6F4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7F1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7F2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7F3.CpG_merged.cov")

file.list.dev.mal = list("../../../mergedCov/dev/DAC3M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC3M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5M3.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC5M4.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC6M3.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DAC7M3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC2M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC2M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC3M4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC4M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC5M3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC6M4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7M1.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7M2.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7M3.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7M4.CpG_merged.cov",
                     "../../../mergedCov/dev/DC7M5.CpG_merged.cov")

#create tabix file
myobj.all=methRead(file.list.dev.all,
               sample.id=list("DAC2F4","DAC2F5","DAC2F6","DAC3F1","DAC3F2","DAC3M1",
                              "DAC3M2","DAC4F1","DAC4F2","DAC4F3","DAC4F4","DAC4F5",
                              "DAC5F1","DAC5F2","DAC5F4","DAC5F5","DAC5M1","DAC5M2",
                              "DAC5M3","DAC5M4","DAC6F1","DAC6F2","DAC6F3","DAC6M1",
                              "DAC6M2","DAC6M3", "DAC7F1","DAC7F2","DAC7F3","DAC7F4",
                              "DAC7F5","DAC7F6","DAC7M1","DAC7M2","DAC7M3",
                              "DC2F1","DC2F2", "DC2M1","DC2M2","DC3F1","DC3F2","DC3F3",
                              "DC3M1","DC3M2","DC3M4","DC4F1","DC4F2","DC4F3","DC4F4",
                              "DC4F5","DC4M1","DC4M2","DC5F1","DC5F2","DC5F3","DC5F4",
                              "DC5F5","DC5M1","DC5M2","DC5M3","DC6F1","DC6F2","DC6F3",
                              "DC6F4","DC6M1","DC6M2","DC6M4","DC7F1","DC7F2","DC7F3",
                              "DC7M1","DC7M2","DC7M3","DC7M4","DC7M5"),
               assembly="guppyWGBS_dev_final",
               pipeline="bismarkCoverage",
               treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
               context="CpG",
              # dbtype = "tabix",
               dbdir = "guppy_dev_DB",
               mincov = 1
)

myobj.fem=methRead(file.list.dev.fem,
               sample.id=list("DAC2F4","DAC2F5","DAC2F6","DAC3F1","DAC3F2",
                              "DAC4F1","DAC4F2","DAC4F3","DAC4F4","DAC4F5",
                              "DAC5F1","DAC5F2","DAC5F4","DAC5F5",
                              "DAC6F1","DAC6F2","DAC6F3",
                              "DAC7F1","DAC7F2","DAC7F3","DAC7F4",
                              "DAC7F5","DAC7F6","DC2F1","DC2F2","DC3F1","DC3F2","DC3F3",
                              "DC4F1","DC4F2","DC4F3","DC4F4",
                              "DC4F5","DC5F1","DC5F2","DC5F3","DC5F4",
                              "DC5F5","DC6F1","DC6F2","DC6F3",
                              "DC6F4","DC7F1","DC7F2","DC7F3"),
               assembly="guppyWGBS_dev_fem",
               pipeline="bismarkCoverage",
               treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
               context="CpG",
              # dbtype = "tabix",
               dbdir = "guppy_dev_DB_fem",
               mincov = 1
)

myobj.mal=methRead(file.list.dev.mal,
               sample.id=list("DAC3M1","DAC3M2","DAC5M1","DAC5M2",
                              "DAC5M3","DAC5M4","DAC6M1",
                              "DAC6M2","DAC6M3", "DAC7M1","DAC7M2","DAC7M3",
                              "DC2M1","DC2M2","DC3M1","DC3M2","DC3M4",
                              "DC4M1","DC4M2","DC5M1","DC5M2","DC5M3",
                              "DC6M1","DC6M2","DC6M4",
                              "DC7M1","DC7M2","DC7M3","DC7M4","DC7M5"),
               assembly="guppyWGBS_dev_mal",
               pipeline="bismarkCoverage",
               treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
               context="CpG",
              # dbtype = "tabix",
               dbdir = "guppy_dev_DB_mal",
               mincov = 1
)

#get coverage stats 
getCoverageStats(myobj.all[[2]], both.strands = FALSE)
getCoverageStats(myobj.fem[[2]], both.strands = FALSE)
getCoverageStats(myobj.mal[[2]], both.strands = FALSE)

#filter out sites in the 99.9th percentile of coverage (PCR bias) and 3x coverage
myobj.all.3X=filterByCoverage(myobj.all,lo.count=3,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9, suffix = "5X_merged")

myobj.fem.3X=filterByCoverage(myobj.fem,lo.count=3,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9, suffix = "5X_merged")

myobj.mal.3X=filterByCoverage(myobj.mal,lo.count=3,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9, suffix = "5X_merged")

#normalize by median coverage
norm.myobj.all.3X=normalizeCoverage(myobj.all.3X, method="median")
norm.myobj.fem.3X=normalizeCoverage(myobj.fem.3X, method="median")
norm.myobj.mal.3X=normalizeCoverage(myobj.mal.3X, method="median")

## Remove sex chromosomes and unplacex scaffolds ##
#prepare GRanges object for chromosomes to keep 
keep.chr.noXY <- GRanges(seqnames = c("NC_024331.1", "NC_024332.1", "NC_024333.1", "NC_024334.1", "NC_024335.1",
                                      "NC_024336.1", "NC_024337.1", "NC_024338.1", "NC_024339.1", "NC_024340.1",
                                      "NC_024341.1", "NC_024343.1", "NC_024344.1", "NC_024345.1", "NC_024346.1",
                                      "NC_024347.1", "NC_024348.1", "NC_024349.1", "NC_024350.1", "NC_024351.1",
                                      "NC_024352.1", "NC_024353.1"),
                         ranges=IRanges(start = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                        end = c(34115677,46286544,35265442,31497199,33908744,
                                                31529174,31413364,27946405,34117797,32819797,
                                                28875558,33524197,28338960,30644713,33199742,
                                                30788009,22026651,28470737,26385442,25773841,
                                                25248790,18084596)),
                         strand="*")

seqlengths(keep.chr.noXY)=c(34115677,46286544,35265442,31497199,33908744,
                            31529174,31413364,27946405,34117797,32819797,
                            28875558,33524197,28338960,30644713,33199742,
                            30788009,22026651,28470737,26385442,25773841,
                            25248790,18084596)

#to remove unplaced scaffolds
keep.chr.allchr <- GRanges(seqnames = c("NC_024331.1", "NC_024332.1", "NC_024333.1", "NC_024334.1", "NC_024335.1",
                                        "NC_024336.1", "NC_024337.1", "NC_024338.1", "NC_024339.1", "NC_024340.1",
                                        "NC_024341.1", "NC_024342.1", "NC_024343.1", "NC_024344.1", "NC_024345.1", "NC_024346.1",
                                        "NC_024347.1", "NC_024348.1", "NC_024349.1", "NC_024350.1", "NC_024351.1",
                                        "NC_024352.1", "NC_024353.1"),
                           ranges=IRanges(start = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                          end = c(34115677,46286544,35265442,31497199,33908744,
                                                  31529174,31413364,27946405,34117797,32819797,
                                                  28875558,26439574,33524197,28338960,30644713,33199742,
                                                  30788009,22026651,28470737,26385442,25773841,
                                                  25248790,18084596)),
                           strand="*")

seqlengths(keep.chr.allchr)=c(34115677,46286544,35265442,31497199,33908744,
                              31529174,31413364,27946405,34117797,32819797,
                              28875558,26439574,33524197,28338960,30644713,33199742,
                              30788009,22026651,28470737,26385442,25773841,
                              25248790,18084596)


myobj3X.all.subset <- selectByOverlap(norm.myobj.all.3X, keep.chr.noXY)
myobj3X.fem.subset <- selectByOverlap(norm.myobj.fem.3X, keep.chr.allchr)
myobj3X.mal.subset <- selectByOverlap(norm.myobj.mal.3X, keep.chr.allchr)

## Find DMRs ##
#read in tank covariate data
covariates.all <- data.frame(tank=c("AC2","AC2","AC2","AC3","AC3","AC3",
                                "AC3","AC4","AC4","AC4","AC4","AC4",
                                "AC5","AC5","AC5","AC5","AC5",
                                "AC5","AC5","AC5","AC6","AC6","AC6",
                                "AC6","AC6","AC6","AC7","AC7","AC7",
                                "AC7","AC7","AC7","AC7","AC7","AC7",
                                "C2","C2","C2","C2","C3","C3",
                                "C3","C3","C3","C3","C4","C4",
                                "C4","C4","C4","C4","C4","C5",
                                "C5","C5","C5","C5","C5","C5",
                                "C5","C6","C6","C6","C6","C6",
                                "C6","C6","C7","C7","C7","C7",
                                "C7","C7","C7","C7"),
                        sex=c("F","F","F","F","F","M",
                              "M","F","F","F","F","F",
                              "F","F","F","F","M","M",
                              "M","M","F","F","F","M",
                              "M","M", "F","F","F","F",
                              "F","F","M","M","M",
                              "F","F", "M","M","F","F","F",
                              "M","M","M","F","F","F","F",
                              "F","M","M","F","F","F","F",
                              "F","M","M","M","F","F","F",
                              "F","M","M","M","F","F","F",
                              "M","M","M","M","M"),
                         stringsAsFactors = TRUE)

covariates.fem <- data.frame(tank=c("AC2","AC2","AC2","AC3","AC3",
                                "AC4","AC4","AC4","AC4","AC4",
                                "AC5","AC5","AC5","AC5",
                                "AC6","AC6","AC6",
                                "AC7","AC7","AC7",
                                "AC7","AC7","AC7",
                                "C2","C2","C3","C3",
                                "C3","C4","C4",
                                "C4","C4","C4","C5",
                                "C5","C5","C5","C5",
                                "C6","C6","C6","C6",
                                "C7","C7","C7"),
                         stringsAsFactors = TRUE)

covariates.mal <- data.frame(tank=c("AC3","AC3",
                                "AC5","AC5","AC5","AC5",
                                "AC6","AC6","AC6",
                                "AC7","AC7","AC7",
                                "C2","C2","C3","C3",
                                "C3","C4","C4",
                                "C5","C5","C5",
                                "C6","C6","C6",
                                "C7","C7","C7","C7",
                                "C7"),
                         stringsAsFactors = TRUE)

#tile into 100 bp windows with min coverage 10X and 5X
tiles.5X.all <- tileMethylCounts(myobj3X.all.subset, win.size = 100, step.size = 100, cov.bases = 5)
tiles.5X.fem <- tileMethylCounts(myobj3X.fem.subset, win.size = 100, step.size = 100, cov.bases = 5)
tiles.5X.mal <- tileMethylCounts(myobj3X.mal.subset, win.size = 100, step.size = 100, cov.bases = 5)

#check number of tiles 
tiles.5X.all
tiles.5X.fem
tiles.5X.mal

#unite calls for 60% of samples
DMRmeth5X.all <- unite(tiles.5X.all, min.per.group=21L, save.db = FALSE)
DMRmeth5X.fem <- unite(tiles.5X.fem, min.per.group=13L, save.db = FALSE)
DMRmeth5X.mal <- unite(tiles.5X.mal, min.per.group=7L, save.db = FALSE)

#check number of regions retained 
DMRmeth5X.all
DMRmeth5X.fem
DMRmeth5X.mal

#calculate differential methylation 
DMRmyDiff5X.all <- calculateDiffMeth(DMRmeth5X.all, mc.cores=2, test="Chisq", covariates=covariates.all, save.db = TRUE, suffix = "myDiffDMR5X_all")
DMRmyDiff5X.fem <- calculateDiffMeth(DMRmeth5X.fem, mc.cores=2, test="Chisq", covariates=covariates.fem, save.db = TRUE, suffix = "myDiffDMR5X_fem")
DMRmyDiff5X.mal <- calculateDiffMeth(DMRmeth5X.mal, mc.cores=2, test="Chisq", covariates=covariates.mal, save.db = TRUE, suffix = "myDiffDMR5X_mal")

#call significant methylation
DMRdiffMeth5X.all <- getMethylDiff(DMRmyDiff5X.all, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "diffMethDMR5X_all")
DMRdiffMeth5X.fem <- getMethylDiff(DMRmyDiff5X.fem, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "diffMethDMR5X_fem")
DMRdiffMeth5X.mal <- getMethylDiff(DMRmyDiff5X.mal, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "diffMethDMR5X_mal")

#check number of DMRs
DMRdiffMeth5X.all
DMRdiffMeth5X.fem
DMRdiffMeth5X.mal

#get meth per chromosome
DMRdiffMethChr5X.all <- diffMethPerChr(DMRmyDiff5X.all, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=20, save.db =  TRUE, suffix = "chrDMR5X_all")
DMRdiffMethChr5X.fem <- diffMethPerChr(DMRmyDiff5X.fem, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=20, save.db =  TRUE, suffix = "chrDMR5X_fem")
DMRdiffMethChr5X.mal <- diffMethPerChr(DMRmyDiff5X.mal, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=20, save.db =  TRUE, suffix = "chrDMR5X_mal")

DMRdiffMethChr5X.all
DMRdiffMethChr5X.fem
DMRdiffMethChr5X.mal

## Save R objects ##
saveRDS(tiles.5X.all, file = "./DMR_res/DMRtiles5X_all.RDS")
saveRDS(DMRmeth5X.all, file = "./DMR_res/DMRmeth5X_all.RDS")
saveRDS(DMRmyDiff5X.all, file = "./DMR_res/DMRmydiff5X_all.RDS")
saveRDS(DMRdiffMeth5X.all, file = "./DMR_res/DMRdiffmeth5X_all.RDS")
saveRDS(DMRdiffMethChr5X.all, file = "./DMR_res/DMRdiffmethchr5X_all.RDS")
saveRDS(getData(DMRmeth5X.all), file = "./DMR_res/DMRmeth5X_data_all.RDS")
saveRDS(getData(DMRmyDiff5X.all), file = "./DMR_res/DMRmydiff5X_data_all.RDS")
saveRDS(getData(DMRdiffMeth5X.all), file = "./DMR_res/DMRdiffmeth5X_data_all.RDS")

saveRDS(tiles.5X.fem, file = "./DMR_res/DMRtiles5X_fem.RDS")
saveRDS(DMRmeth5X.fem, file = "./DMR_res/DMRmeth5X_fem.RDS")
saveRDS(DMRmyDiff5X.fem, file = "./DMR_res/DMRmydiff5X_fem.RDS")
saveRDS(DMRdiffMeth5X.fem, file = "./DMR_res/DMRdiffmeth5X_fem.RDS")
saveRDS(DMRdiffMethChr5X.fem, file = "./DMR_res/DMRdiffmethchr5X_fem.RDS")
saveRDS(getData(DMRmeth5X.fem), file = "./DMR_res/DMRmeth5X_data_fem.RDS")
saveRDS(getData(DMRmyDiff5X.fem), file = "./DMR_res/DMRmydiff5X_data_fem.RDS")
saveRDS(getData(DMRdiffMeth5X.fem), file = "./DMR_res/DMRdiffmeth5X_data_fem.RDS")

saveRDS(tiles.5X.mal, file = "./DMR_res/DMRtiles5X_mal.RDS")
saveRDS(DMRmeth5X.mal, file = "./DMR_res/DMRmeth5X_mal.RDS")
saveRDS(DMRmyDiff5X.mal, file = "./DMR_res/DMRmydiff5X_mal.RDS")
saveRDS(DMRdiffMeth5X.mal, file = "./DMR_res/DMRdiffmeth5X_mal.RDS")
saveRDS(DMRdiffMethChr5X.mal, file = "./DMR_res/DMRdiffmethchr5X_mal.RDS")
saveRDS(getData(DMRmeth5X.mal), file = "./DMR_res/DMRmeth5X_data_mal.RDS")
saveRDS(getData(DMRmyDiff5X.mal), file = "./DMR_res/DMRmydiff5X_data_mal.RDS")
saveRDS(getData(DMRdiffMeth5X.mal), file = "./DMR_res/DMRdiffmeth5X_data_mal.RDS")







