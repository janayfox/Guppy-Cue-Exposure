###########################################################################################
### Goal: Remove sex chr and unplaced scaffolds and methylation calling for DMRs 
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

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X and 5X
tiles.10X <- tileMethylCounts(myobj.3X, win.size = 100, step.size = 100, cov.bases = 10, save.db = TRUE, suffix = "tiles10X")
tiles.5X <- tileMethylCounts(myobj.3X, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles5X")

#remove sex chr (LG12) and unplaced scaffolds
tiles.10X.subset <- selectByOverlap(tiles.10X, keep.chr)
tiles.5X.subset <- selectByOverlap(tiles.5X, keep.chr)

#convert back to DB
tiles.10X.subsetDB <- makeMethylDB(tiles.10X.subset, "tiles_subset_10XDB")
tiles.5X.subsetDB <- makeMethylDB(tiles.5X.subset, "tiles_subset_5XDB")

#unite calls
DMR.meth.2L10X <- unite(tiles.10X.subsetDB, min.per.group=2L, save.db = TRUE, suffix = "DMR_unite_2L10X")
DMR.meth.10L10X <- unite(tiles.10X.subsetDB, min.per.group=10L, save.db = TRUE, suffix = "DMR_unite_10L10X")
DMR.meth.2L5X <- unite(tiles.5X.subsetDB, min.per.group=2L, save.db = TRUE, suffix = "DMR_unite_2L5X")
DMR.meth.10L5X <- unite(tiles.5X.subsetDB, min.per.group=10L, save.db = TRUE, suffix = "DMR_unite_10L5X")

#calculate differential methylation 
DMR.myDiff.2L10X <- calculateDiffMeth(DMR.meth.2L10X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "DMR_myDiff_2L10X")
DMR.myDiff.10L10X <- calculateDiffMeth(DMR.meth.10L10X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "DMR_myDiff_10L10X")
DMR.myDiff.2L5X <- calculateDiffMeth(DMR.meth.2L5X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "DMR_myDiff_2L5X")
DMR.myDiff.10L5X <- calculateDiffMeth(DMR.meth.10L5X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "DMR_myDiff_10L5X")
 
#call significant methylation
DMR.diffMeth.2L10X <- getMethylDiff(DMR.myDiff.2L10X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_2L10X")
DMR.diffMeth.10L10X <- getMethylDiff(DMR.myDiff.10L10X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_10L10X")
DMR.diffMeth.2L5X <- getMethylDiff(DMR.myDiff.2L5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_2L5X")
DMR.diffMeth.10L5X <- getMethylDiff(DMR.myDiff.10L5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_10L5X")

#get meth per chromosome
DMR.diffMethChr.2L10X <- diffMethPerChr(DMR.myDiff.2L10X, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "DMR_chr_2L10X")
DMR.diffMethChr.10L10X <- diffMethPerChr(DMR.myDiff.10L10X, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "DMR_chr_10L10X")
DMR.diffMethChr.2L5X <- diffMethPerChr(DMR.myDiff.2L5X, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "DMR_chr_2L5X")
DMR.diffMethChr.10L5X <- diffMethPerChr(DMR.myDiff.10L5X, plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "DMR_chr_10L5X")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/06_methylCallingDMR-backup.RData")

q(save="yes")


