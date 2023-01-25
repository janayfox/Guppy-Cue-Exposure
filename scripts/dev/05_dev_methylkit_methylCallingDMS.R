###########################################################################################
### Goal: Remove sex chr and unplaced scaffolds and methylation calling for DMSs 
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

setwd("/scratch/janayfox/guppyWGBS/")
load(file=".RData")

#prepare GRanges object for chromosomes to keep
keep.chr <- GRanges(seqnames = c("NC_024331.1", "NC_024332.1", "NC_024333.1", "NC_024334.1", "NC_024335.1",
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

seqlengths(keep.chr)=c(34115677,46286544,35265442,31497199,33908744,
                       31529174,31413364,27946405,34117797,32819797,
                       28875558,33524197,28338960,30644713,33199742,
                       30788009,22026651,28470737,26385442,25773841,
                       25248790,18084596)

#read in tank covariate data
covariates <- data.frame(tank=c("AC2","AC2","AC2","AC3","AC3","AC3",
                                "AC3","AC4","AC4","AC4","AC4","AC4","AC4",
                                "AC5", "AC5","AC5","AC5","AC5","AC5","AC5",
                                "AC5","AC6","AC6","AC6","AC6","AC6", "AC6",
                                "AC7","AC7","AC7","AC7","AC7","AC7","AC7",
                                "AC7","AC7","C2","C2", "C2","C2","C3","C3","C3",
                                "C3","C3","C3","C4","C4","C4","C4",
                                "C4","C4","C4","C5","C5","C5","C5",
                                "C5","C5","C5","C5","C6","C6","C6",
                                "C6","C6","C6","C6","C7","C7","C7",
                                "C7","C7","C7","C7","C7"), stringsAsFactors = TRUE)

 ## Find DMS ##
#remove sex chr (LG12) and unplaced scaffolds
myobj.3X.subset <- selectByOverlap(myobj.3X, keep.chr)

#convert back to DB
myobj.3X.subsetDB <- makeMethylDB(myobj.3X.subset, "subset_3X_DB")

#filter out sites for 10x coverage for DMS
myobj.10X.subset=filterByCoverage(myobj.3X.subsetDB,lo.count=10,lo.perc=NULL,
                            hi.count=NULL, hi.perc=NULL, suffix = "subset_10X")

#also try 5X coverage
myobj.5X.subset=filterByCoverage(myobj.3X.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "subset_5X")

#normalize by median coverage
norm.myobj.10X.subset=normalizeCoverage(myobj.10X.subset, method="median")
norm.myobj.5X.subset=normalizeCoverage(myobj.5X.subset, method="median")

#keep sites for min of 2 and 10 fish per group
meth.2L.10X=unite(norm.myobj.10X.subset, min.per.group=2L, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_2L10X")
meth.10L.10X=unite(norm.myobj.10X.subset, min.per.group=10L, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_10L10X")
meth.2L.5X=unite(norm.myobj.5X.subset, min.per.group=2L, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_2L5X")
meth.10L.5X=unite(norm.myobj.5X.subset, min.per.group=10L, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_10L5X")

# calculate differential methylation
DMS.myDiff.2L.10X <- calculateDiffMeth(meth.2L.10X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "DMS_myDiff_2L10X")
DMS.myDiff.10L.10X <- calculateDiffMeth(meth.10L.10X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "DMS_myDiff_10L10X")
DMS.myDiff.2L.5X <- calculateDiffMeth(meth.2L.5X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "DMS_myDiff_2L5X")
DMS.myDiff.10L.5X <- calculateDiffMeth(meth.10L.5X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "DMS_myDiff_10L5X")

#call significant methylation
DMS.diffMeth.2L.10X <- getMethylDiff(DMS.myDiff.2L.10X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_2L10X")
DMS.diffMeth.10L.10X <- getMethylDiff(DMS.myDiff.10L.10X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_10L10X")
DMS.diffMeth.2L.5X <- getMethylDiff(DMS.myDiff.2L.5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_2L5X")
DMS.diffMeth.10L.5X <- getMethylDiff(DMS.myDiff.10L.5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_10L5X")

#get meth per chromosome
DMS.diffMethChr.2L10X <- diffMethPerChr(DMS.myDiff.2L.10X,plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "DMS_chr_2L10X")
DMS.diffMethChr.10L10X <- diffMethPerChr(DMS.myDiff.10L.10X,plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "DMS_chr_10L10X")
DMS.diffMethChr.2L5X <- diffMethPerChr(DMS.myDiff.2L.5X,plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "DMS_chr_2L5X")
DMS.diffMethChr.10L5X <- diffMethPerChr(DMS.myDiff.10L.5X,plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "DMS_chr_10L5X")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/05_methylCallingDMS-backup.RData")

q(save="yes")
