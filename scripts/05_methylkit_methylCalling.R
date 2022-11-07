###########################################################################################
### Goal: Remove sex chr and unplaced scaffolds and methylation calling for DMRs and DMSs 
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

#prepare GRanges object
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
covariates=data.frame(tank=c("AC2","AC2","AC2","AC3","AC3","AC3",
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
#subset
#norm.10X.subset <- selectByOverlap(norm.myobj.10X[[1]], keep.chr)

#convert back to DB
#norm.10X.subsetDB <- makeMethylDB(norm.10X.subset, "10XsubsetDB")

#keep sites for min of 10 fish per group
#meth10X=unite(norm10X.subsetDB, min.per.group=10L, destrand=TRUE, save.db = TRUE, suffix = "DMS_10L_unite_subset")
#doesnt work.. but it does work on norm.myobj.10x -> it has become a methylRaw not rawlist .. not sure how to fix this 

#calculate diff meth
#myDiff=calculateDiffMeth(meth10X, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "DMS_10L_myDiff")
#diffMeth <- getMethylDiff(myDiff, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_10L_diffMeth")
#diffMethChr <- diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "DMS_10L_Chr")

#myDiff_df <- getData(myDiff)

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X
tiles <- tileMethylCounts(norm.myobj.3X, win.size = 100, step.size = 100, cov.bases = 10, save.db = TRUE, suffix = "tiles_100")
#remove sex chr (LG12) and unplaced scaffolds
tiles.subset <- selectByOverlap(tiles, keep.chr)
tiles.subsetDB <- makeMethylDB(tiles.subset, "10XsubsetDB")
#unite calls
DMR.meth <- unite(tiles.subsetDB, min.per.group=10L, destrand=TRUE, save.db = TRUE, suffix = "DMR_10L_unite")
#calculate differential methylation 
DMR.myDiff=calculateDiffMeth(DMR.meth, mc.cores=1, covariates=covariates, save.db = TRUE, suffix = "myDiff")
#call significant methylation
DMR.diffMeth <- getMethylDiff(DMR.myDiff, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "diffMeth")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/05_methylCalling-backup.RData")

q(save="yes")


