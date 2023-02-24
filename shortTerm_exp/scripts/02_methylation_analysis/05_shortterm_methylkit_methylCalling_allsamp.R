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

setwd("/scratch/janayfox/guppyWGBS_shortterm/all_samp/")
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
covariates.all <- data.frame(tank=c("AC10","AC10","AC11","AC11","AC13",
                                    "AC13","AC14","AC14","AC15","AC15",
                                    "AC16","AC16","AC1","AC1","AC2",
                                    "AC2","AC3","AC3","AC4","AC4",
                                    "AC5","AC5","AC6","AC6","AC7",
                                    "AC7","AC8","AC8","AC9","AC9",
                                    "C10","C10","C11","C11","C12",
                                    "C12","C13","C13","C14","C14",
                                    "C15","C15","C1","C1","C2",
                                    "C2","C3","C3","C4","C4",
                                    "C5","C5","C6","C6","C7",
                                    "C7","C8","C8","C9", "C9"), 
                             sex=c("F","M","F","M","F","M","F","M","F","M",
                                   "F","M","F","M","F","M","F","M","F","M",
                                   "F","M","F","M","F","M","F","M","F","M",
                                   "F","M","F","M","F","M","F","M","F","M",
                                   "F","M","F","M","F","M","F","M","F","M",
                                   "F","M","F","M","F","M","F","M","F", "M"), 
                             stringsAsFactors = TRUE)

## Find DMS ##
#remove sex chr (LG12) and unplaced scaffolds
myobj_all.subset <- selectByOverlap(myobj_all.filt, keep.chr.noXY)
myobj_all.fem.subset <- selectByOverlap(myobj_all.fem.filt, keep.chr.allchr)
myobj_all.mal.subset <- selectByOverlap(myobj_all.mal.filt, keep.chr.allchr)

#convert back to DB
myobj_all.subsetDB <- makeMethylDB(myobj_all.subset, "all_subset_DB")
myobj_all.fem.subsetDB <- makeMethylDB(myobj_all.fem.subset, "all_fem_subset_DB")
myobj_all.mal.subsetDB <- makeMethylDB(myobj_all.mal.subset, "all_mal_subset_DB")

#filter out sites for 5x coverage 
myobj_all.5X=filterByCoverage(myobj_all.subsetDB,lo.count=5,lo.perc=NULL,
                              hi.count=NULL, hi.perc=NULL, suffix = "all_5X")
myobj_all.fem.5X=filterByCoverage(myobj_all.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "all_fem_5X")
myobj_all.mal.5X=filterByCoverage(myobj_all.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "all_mal_5X")

#normalize by median coverage
myobj_all.norm=normalizeCoverage(myobj_all.5X, method="all_median")
myobj_all.fem.norm=normalizeCoverage(myobj_all.fem.5X, method="all_fem_median")
myobj_all.mal.norm=normalizeCoverage(myobj_all.mal.5X, method="all_mal_median")

#unite sites 
meth.all=unite(myobj_all.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_all")
meth.all.fem=unite(myobj_all.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_all_fem")
meth.all.mal=unite(myobj_all.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_all_mal")

#calculate differential methylation
DMS.myDiff.all <- calculateDiffMeth(meth.all, covariates=covariates.all, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_all")
DMS.myDiff.all.fem <- calculateDiffMeth(meth.all.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_all_fem")
DMS.myDiff.all.mal <- calculateDiffMeth(meth.all.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_all_mal")

#call significant methylation
DMS.diffMeth.all <- getMethylDiff(DMS.myDiff.all, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_all")
DMS.diffMeth.all.fem <- getMethylDiff(DMS.myDiff.all.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_all.fem")
DMS.diffMeth.all.mal <- getMethylDiff(DMS.myDiff.all.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_all.mal")

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X and 5X
tiles.all <- tileMethylCounts(myobj_all.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_all")
tiles.all.fem <- tileMethylCounts(myobj_all.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_allfem")
tiles.all.mal <- tileMethylCounts(myobj_all.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_allmal")

#unite calls
DMR.meth.all <- unite(tiles.all, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_all")
DMR.meth.all.fem <- unite(tiles.all.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_allfem")
DMR.meth.all.mal <- unite(tiles.all.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_allmal")

#calculate differential methylation 
DMR.myDiff.all <- calculateDiffMeth(DMR.meth.all, mc.cores=1, covariates=covariates.all, save.db = TRUE, suffix = "DMR_myDiff_all")
DMR.myDiff.all.fem <- calculateDiffMeth(DMR.meth.all.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_all")
DMR.myDiff.all.mal <- calculateDiffMeth(DMR.meth.all.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_all")

#call significant methylation
DMR.diffMeth.all <- getMethylDiff(DMR.myDiff.all, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_all")
DMR.diffMeth.all.fem <- getMethylDiff(DMR.myDiff.all.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_allfem")
DMR.diffMeth.all.mal <- getMethylDiff(DMR.myDiff.all.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_allmal")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/05_shortterm_methylCalling_allsamp-backup.RData")

q(save="yes")

