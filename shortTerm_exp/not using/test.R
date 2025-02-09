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

setwd("/scratch/janayfox/guppyWGBS_shortterm/1h/")

## Prepare tabix files
#create lists of file locations
file.list.1h = list("../mergedCov/st/ST2AC14F.CpG_merged.cov",
                    "../mergedCov/st/ST2AC14M.CpG_merged.cov",
                    "../mergedCov/st/ST2AC2F.CpG_merged.cov",
                    "../mergedCov/st/ST2AC2M.CpG_merged.cov",
                    "../mergedCov/st/ST2AC9F.CpG_merged.cov",
                    "../mergedCov/st/ST2AC9M.CpG_merged.cov",
                    "../mergedCov/st/ST2C14F.CpG_merged.cov",
                    "../mergedCov/st/ST2C14M.CpG_merged.cov",
                    "../mergedCov/st/ST2C3F.CpG_merged.cov",
                    "../mergedCov/st/ST2C3M.CpG_merged.cov",
                    "../mergedCov/st/ST2C9F.CpG_merged.cov",
                    "../mergedCov/st/ST2C9M.CpG_merged.cov")

file.list.1h.v2 = list("../results/bismark_methylation_calls/stranded_CpG_report/ST2AC14F.CpG_report.txt.gz",
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
                 

#create tabix file
myobj.1h=methRead(file.list.1h,
                  sample.id=list("ST2AC14F","ST2AC14M","ST2AC2F","ST2AC2M","ST2AC9F","ST2AC9M",
                                 "ST2C14F","ST2C14M","ST2C3F","ST2C3M","ST2C9F", "ST2C9M"),
                  assembly="guppyWGBS_shortterm",
                  pipeline="bismarkCoverage",
                  treatment=c(1,1,1,1,1,1,
                              0,0,0,0,0,0),
                  context="CpG",
                  mincov = 1,
                  dbtype = "tabix",
                  dbdir = "shortterm_1h_DB_merged"
)

myobj.1h.v2=methRead(file.list.1h.v2,
                  sample.id=list("ST2AC14F","ST2AC14M","ST2AC2F","ST2AC2M","ST2AC9F","ST2AC9M",
                                 "ST2C14F","ST2C14M","ST2C3F","ST2C3M","ST2C9F", "ST2C9M"),
                  assembly="guppyWGBS_shortterm",
                  pipeline="bismarkCytosineReport",
                  treatment=c(1,1,1,1,1,1,
                              0,0,0,0,0,0),
                  context="CpG",
                  mincov = 1,
                  dbtype = "tabix",
                  dbdir = "shortterm_1h_DB_not_merged"
)

#get coverage stats 
getCoverageStats(myobj.1h[[2]], both.strands = FALSE)
getCoverageStats(myobj.1h.v2[[2]], both.strands = TRUE)

#filter out sites in the 99.9th percentile of coverage (PCR bias) 
myobj.1h.5X=filterByCoverage(myobj.1h,lo.count=5,lo.perc=NULL,
                               hi.count=NULL, hi.perc=99.9, suffix = "5X")

myobj.1h.5X.v2=filterByCoverage(myobj.1h.v2,lo.count=5,lo.perc=NULL,
                               hi.count=NULL, hi.perc=99.9, suffix = "5X")

#normalize by median coverage
norm.myob.1h.5X=normalizeCoverage(myobj.1h.5X, method="median")
norm.myob.1h.5X.v2=normalizeCoverage(myobj.1h.5X.v2, method="median")

#prepare GRanges object for chromosomes to keep 
#to remove unplaced scaffolds and sex chromosomes
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

myobj.1h.subset <- selectByOverlap(norm.myobj.1h.5X, keep.chr.noXY)
myobj.1h.subset.v2 <- selectByOverlap(norm.myobj.1h.5X.v2, keep.chr.noXY)

##Find DMS##
#unite sites 
DMS.meth.1h.5X=unite(myobj.1h.subset, min.per.group = 2L, destrand=FALSE, save.db = TRUE, suffix = "DMS_unite_1h")
DMS.meth.1h.5X.v2=unite(myobj.1h.subset.v2, min.per.group = 2L, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_1h_v2")

# Check number of CpGs 
DMSmeth.1h.5X
DMSmeth.1h.5X.v2

#enter covariates 
covariates.1h <- data.frame(tank=c("AC14","AC14","AC2","AC2","AC9","AC9",
                                   "C14","C14","C3","C3","C9","C9"), 
                            sex=c("F","M","F","M","F","M",
                                  "F","M","F","M","F","M"), 
                            stringsAsFactors = TRUE)

#calculate differential methylation
DMS.myDiff.1h.5X <- calculateDiffMeth(DMS.meth.1h.5X, covariates=covariates.1h, mc.cores=2, overdispersion="MN", test="Chisq", save.db = TRUE, suffix = "myDiff")
DMS.myDiff.1h.5X.v2 <- calculateDiffMeth(DMS.meth.1h.5X.v2, covariates=covariates.1h, mc.cores=2, overdispersion="MN", test="Chisq", save.db = TRUE, suffix = "myDiff")

#call significant methylation
DMS.diffMeth.1h.5X <- getMethylDiff(DMS.myDiff.1h.5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "diffMeth")
DMS.diffMeth.1h.5X.v2 <- getMethylDiff(DMS.myDiff.1h.5X.v2, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "diffMeth")

#check number of significant DMS
DMS.diffMeth.1h.5X
DMS.diffMeth.1h.5X.v2


