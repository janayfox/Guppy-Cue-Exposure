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

setwd("/scratch/janayfox/guppyWGBS_shortterm/05hr/")
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
covariates.05h <- data.frame(tank=c("AC10","AC10","AC16","AC16","AC3","AC3",
                                    "C10","C10","C13","C13","C3","C3"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

## Find DMS ##
#remove sex chr (LG12) and unplaced scaffolds
myobj_05h.subset <- selectByOverlap(myobj_05h.filt, keep.chr.noXY)
myobj_05h.fem.subset <- selectByOverlap(myobj_05h.fem.filt, keep.chr.allchr)
myobj_05h.mal.subset <- selectByOverlap(myobj_05h.mal.filt, keep.chr.allchr)

#convert back to DB
myobj_05h.subsetDB <- makeMethylDB(myobj_05h.subset, "05h_subset_DB")
myobj_05h.fem.subsetDB <- makeMethylDB(myobj_05h.fem.subset, "05h_fem_subset_DB")
myobj_05h.mal.subsetDB <- makeMethylDB(myobj_05h.mal.subset, "05h_mal_subset_DB")

#filter out sites for 5x coverage 
myobj_05h.5X=filterByCoverage(myobj_05h.subsetDB,lo.count=5,lo.perc=NULL,
                              hi.count=NULL, hi.perc=NULL, suffix = "05h_5X")
myobj_05h.fem.5X=filterByCoverage(myobj_05h.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "05h_fem_5X")
myobj_05h.mal.5X=filterByCoverage(myobj_05h.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "05h_mal_5X")

#normalize by median coverage
myobj_05h.norm=normalizeCoverage(myobj_05h.5X, method="median")
myobj_05h.fem.norm=normalizeCoverage(myobj_05h.fem.5X, method="median")
myobj_05h.mal.norm=normalizeCoverage(myobj_05h.mal.5X, method="median")

#unite sites 
meth.05h=unite(myobj_05h.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_05h")
meth.05h.fem=unite(myobj_05h.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_05h_fem")
meth.05h.mal=unite(myobj_05h.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_05h_mal")

#calculate differential methylation
DMS.myDiff.05h <- calculateDiffMeth(meth.05h, covariates=covariates.05h, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_05h")
DMS.myDiff.05h.fem <- calculateDiffMeth(meth.05h.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_05h_fem")
DMS.myDiff.05h.mal <- calculateDiffMeth(meth.05h.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_05h_mal")

#call significant methylation
DMS.diffMeth.05h <- getMethylDiff(DMS.myDiff.05h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_05h")
DMS.diffMeth.05h.fem <- getMethylDiff(DMS.myDiff.05h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_05h.fem")
DMS.diffMeth.05h.mal <- getMethylDiff(DMS.myDiff.05h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_05h.mal")

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X and 5X
tiles.05h <- tileMethylCounts(myobj_05h.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_05h")
tiles.05h.fem <- tileMethylCounts(myobj_05h.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_05hfem")
tiles.05h.mal <- tileMethylCounts(myobj_05h.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_05hmal")

#unite calls
DMR.meth.05h <- unite(tiles.05h, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_05h")
DMR.meth.05h.fem <- unite(tiles.05h.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_05hfem")
DMR.meth.05h.mal <- unite(tiles.05h.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_05hmal")

#calculate differential methylation 
DMR.myDiff.05h <- calculateDiffMeth(DMR.meth.05h, mc.cores=1, covariates=covariates.05h, save.db = TRUE, suffix = "DMR_myDiff_05h")
DMR.myDiff.05h.fem <- calculateDiffMeth(DMR.meth.05h.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_05h")
DMR.myDiff.05h.mal <- calculateDiffMeth(DMR.meth.05h.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_05h")

#call significant methylation
DMR.diffMeth.05h <- getMethylDiff(DMR.myDiff.05h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05h")
DMR.diffMeth.05h.fem <- getMethylDiff(DMR.myDiff.05h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05hfem")
DMR.diffMeth.05h.mal <- getMethylDiff(DMR.myDiff.05h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05hmal")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/05_shortterm_methylCalling_05hr-backup.RData")

q(save="yes")

