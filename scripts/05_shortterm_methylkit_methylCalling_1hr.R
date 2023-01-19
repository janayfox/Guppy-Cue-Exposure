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

setwd("/scratch/janayfox/guppyWGBS_shortterm/1hr/")
load(file=".RData")

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

#enter covariates 
covariates.1h <- data.frame(tank=c("AC14","AC14","AC2","AC2","AC9","AC9",
                                   "C14","C14","C3","C3","C9","C9"), 
                            sex=c("F","M","F","M","F","M",
                                  "F","M","F","M","F","M"), 
                            stringsAsFactors = TRUE)

## Find DMS ##
#remove sex chr (LG12) and unplaced scaffolds
myobj_1h.subset <- selectByOverlap(myobj_1h.filt, keep.chr.noXY)
myobj_1h.fem.subset <- selectByOverlap(myobj_1h.fem.filt, keep.chr.allchr)
myobj_1h.mal.subset <- selectByOverlap(myobj_1h.mal.filt, keep.chr.allchr)

#convert back to DB
myobj_1h.subsetDB <- makeMethylDB(myobj_1h.subset, "1h_subset_DB")
myobj_1h.fem.subsetDB <- makeMethylDB(myobj_1h.fem.subset, "1h_fem_subset_DB")
myobj_1h.mal.subsetDB <- makeMethylDB(myobj_1h.mal.subset, "1h_mal_subset_DB")

#filter out sites for 5x coverage 
myobj_1h.5X=filterByCoverage(myobj_1h.subsetDB,lo.count=5,lo.perc=NULL,
                             hi.count=NULL, hi.perc=NULL, suffix = "1h_5X")
myobj_1h.fem.5X=filterByCoverage(myobj_1h.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                 hi.count=NULL, hi.perc=NULL, suffix = "1h_fem_5X")
myobj_1h.mal.5X=filterByCoverage(myobj_1h.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                 hi.count=NULL, hi.perc=NULL, suffix = "1h_mal_5X")

#normalize by median coverage
myobj_1h.norm=normalizeCoverage(myobj_1h.5X, method="median")
myobj_1h.fem.norm=normalizeCoverage(myobj_1h.fem.5X, method="median")
myobj_1h.mal.norm=normalizeCoverage(myobj_1h.mal.5X, method="median")

#unite sites 
meth.1h=unite(myobj_1h.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_1h")
meth.1h.fem=unite(myobj_1h.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_1h_fem")
meth.1h.mal=unite(myobj_1h.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_1h_mal")

#calculate differential methylation
DMS.myDiff.1h <- calculateDiffMeth(meth.1h, covariates=covariates.1h, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_1h")
DMS.myDiff.1h.fem <- calculateDiffMeth(meth.1h.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_1h_fem")
DMS.myDiff.1h.mal <- calculateDiffMeth(meth.1h.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_1h_mal")

#call significant methylation
DMS.diffMeth.1h <- getMethylDiff(DMS.myDiff.1h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_1h")
DMS.diffMeth.1h.fem <- getMethylDiff(DMS.myDiff.1h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_1h.fem")
DMS.diffMeth.1h.mal <- getMethylDiff(DMS.myDiff.1h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_1h.mal")

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X and 5X
tiles.1h <- tileMethylCounts(myobj_1h.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_1h")
tiles.1h.fem <- tileMethylCounts(myobj_1h.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_1hfem")
tiles.1h.mal <- tileMethylCounts(myobj_1h.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_1hmal")

#unite calls
DMR.meth.1h <- unite(tiles.1h, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_1h")
DMR.meth.1h.fem <- unite(tiles.1h.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_1hfem")
DMR.meth.1h.mal <- unite(tiles.1h.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_1hmal")

#calculate differential methylation 
DMR.myDiff.1h <- calculateDiffMeth(DMR.meth.1h, mc.cores=1, covariates=covariates.1h, save.db = TRUE, suffix = "DMR_myDiff_1h")
DMR.myDiff.1h.fem <- calculateDiffMeth(DMR.meth.1h.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_1h")
DMR.myDiff.1h.mal <- calculateDiffMeth(DMR.meth.1h.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_1h")

#call significant methylation
DMR.diffMeth.1h <- getMethylDiff(DMR.myDiff.1h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_1h")
DMR.diffMeth.1h.fem <- getMethylDiff(DMR.myDiff.1h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_1hfem")
DMR.diffMeth.1h.mal <- getMethylDiff(DMR.myDiff.1h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_1hmal")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/05_shortterm_methylCalling_1hr-backup.RData")

q(save="yes")

