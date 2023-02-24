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

setwd("/scratch/janayfox/guppyWGBS_shortterm/24hr/")
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
covariates.24h <- data.frame(tank=c("AC15","AC15","AC4","AC4","AC8","AC8",
                                    "C15","C15","C4","C4","C8","C8"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

## Find DMS ##
#remove sex chr (LG12) and unplaced scaffolds
myobj_24h.subset <- selectByOverlap(myobj_24h.filt, keep.chr.noXY)
myobj_24h.fem.subset <- selectByOverlap(myobj_24h.fem.filt, keep.chr.allchr)
myobj_24h.mal.subset <- selectByOverlap(myobj_24h.mal.filt, keep.chr.allchr)

#convert back to DB
myobj_24h.subsetDB <- makeMethylDB(myobj_24h.subset, "24h_subset_DB")
myobj_24h.fem.subsetDB <- makeMethylDB(myobj_24h.fem.subset, "24h_fem_subset_DB")
myobj_24h.mal.subsetDB <- makeMethylDB(myobj_24h.mal.subset, "24h_mal_subset_DB")

#filter out sites for 5x coverage 
myobj_24h.5X=filterByCoverage(myobj_24h.subsetDB,lo.count=5,lo.perc=NULL,
                              hi.count=NULL, hi.perc=NULL, suffix = "24h_5X")
myobj_24h.fem.5X=filterByCoverage(myobj_24h.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "24h_fem_5X")
myobj_24h.mal.5X=filterByCoverage(myobj_24h.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "24h_mal_5X")

#normalize by median coverage
myobj_24h.norm=normalizeCoverage(myobj_24h.5X, method="median")
myobj_24h.fem.norm=normalizeCoverage(myobj_24h.fem.5X, method="median")
myobj_24h.mal.norm=normalizeCoverage(myobj_24h.mal.5X, method="median")

#unite sites 
meth.24h=unite(myobj_24h.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_24h")
meth.24h.fem=unite(myobj_24h.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_24h_fem")
meth.24h.mal=unite(myobj_24h.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_24h_mal")

#calculate differential methylation
DMS.myDiff.24h <- calculateDiffMeth(meth.24h, covariates=covariates.24h, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_24h")
DMS.myDiff.24h.fem <- calculateDiffMeth(meth.24h.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_24h_fem")
DMS.myDiff.24h.mal <- calculateDiffMeth(meth.24h.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_24h_mal")

#call significant methylation
DMS.diffMeth.24h <- getMethylDiff(DMS.myDiff.24h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_24h")
DMS.diffMeth.24h.fem <- getMethylDiff(DMS.myDiff.24h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_24h.fem")
DMS.diffMeth.24h.mal <- getMethylDiff(DMS.myDiff.24h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_24h.mal")

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X and 5X
tiles.24h <- tileMethylCounts(myobj_24h.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_24h")
tiles.24h.fem <- tileMethylCounts(myobj_24h.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_24hfem")
tiles.24h.mal <- tileMethylCounts(myobj_24h.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_24hmal")

#unite calls
DMR.meth.24h <- unite(tiles.24h, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_24h")
DMR.meth.24h.fem <- unite(tiles.24h.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_24hfem")
DMR.meth.24h.mal <- unite(tiles.24h.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_24hmal")

#calculate differential methylation 
DMR.myDiff.24h <- calculateDiffMeth(DMR.meth.24h, mc.cores=1, covariates=covariates.24h, save.db = TRUE, suffix = "DMR_myDiff_24h")
DMR.myDiff.24h.fem <- calculateDiffMeth(DMR.meth.24h.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_24h")
DMR.myDiff.24h.mal <- calculateDiffMeth(DMR.meth.24h.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_24h")

#call significant methylation
DMR.diffMeth.24h <- getMethylDiff(DMR.myDiff.24h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_24h")
DMR.diffMeth.24h.fem <- getMethylDiff(DMR.myDiff.24h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_24hfem")
DMR.diffMeth.24h.mal <- getMethylDiff(DMR.myDiff.24h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_24hmal")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/05_shortterm_methylCalling_24hr-backup.RData")

q(save="yes")

