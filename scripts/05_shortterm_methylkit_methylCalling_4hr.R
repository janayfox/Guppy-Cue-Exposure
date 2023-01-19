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

setwd("/scratch/janayfox/guppyWGBS_shortterm/4hr/")
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
covariates.4h <- data.frame(tank=c("AC13","AC13","AC1","AC1","AC7","AC7",
                                   "C12","C12","C1","C1","C7","C7"), 
                            sex=c("F","M","F","M","F","M",
                                  "F","M","F","M","F","M"), 
                            stringsAsFactors = TRUE)

## Find DMS ##
#remove sex chr (LG12) and unplaced scaffolds
myobj_4h.subset <- selectByOverlap(myobj_4h.filt, keep.chr.noXY)
myobj_4h.fem.subset <- selectByOverlap(myobj_4h.fem.filt, keep.chr.allchr)
myobj_4h.mal.subset <- selectByOverlap(myobj_4h.mal.filt, keep.chr.allchr)

#convert back to DB
myobj_4h.subsetDB <- makeMethylDB(myobj_4h.subset, "4h_subset_DB")
myobj_4h.fem.subsetDB <- makeMethylDB(myobj_4h.fem.subset, "4h_fem_subset_DB")
myobj_4h.mal.subsetDB <- makeMethylDB(myobj_4h.mal.subset, "4h_mal_subset_DB")

#filter out sites for 5x coverage 
myobj_4h.5X=filterByCoverage(myobj_4h.subsetDB,lo.count=5,lo.perc=NULL,
                             hi.count=NULL, hi.perc=NULL, suffix = "4h_5X")
myobj_4h.fem.5X=filterByCoverage(myobj_4h.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                 hi.count=NULL, hi.perc=NULL, suffix = "4h_fem_5X")
myobj_4h.mal.5X=filterByCoverage(myobj_4h.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                 hi.count=NULL, hi.perc=NULL, suffix = "4h_mal_5X")

#normalize by median coverage
myobj_4h.norm=normalizeCoverage(myobj_4h.5X, method="median")
myobj_4h.fem.norm=normalizeCoverage(myobj_4h.fem.5X, method="median")
myobj_4h.mal.norm=normalizeCoverage(myobj_4h.mal.5X, method="median")

#unite sites 
meth.4h=unite(myobj_4h.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_4h")
meth.4h.fem=unite(myobj_4h.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_4h_fem")
meth.4h.mal=unite(myobj_4h.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_4h_mal")

#calculate differential methylation
DMS.myDiff.4h <- calculateDiffMeth(meth.4h, covariates=covariates.4h, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_4h")
DMS.myDiff.4h.fem <- calculateDiffMeth(meth.4h.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_4h_fem")
DMS.myDiff.4h.mal <- calculateDiffMeth(meth.4h.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_4h_mal")

#call significant methylation
DMS.diffMeth.4h <- getMethylDiff(DMS.myDiff.4h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_4h")
DMS.diffMeth.4h.fem <- getMethylDiff(DMS.myDiff.4h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_4h.fem")
DMS.diffMeth.4h.mal <- getMethylDiff(DMS.myDiff.4h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_4h.mal")

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X and 5X
tiles.4h <- tileMethylCounts(myobj_4h.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_4h")
tiles.4h.fem <- tileMethylCounts(myobj_4h.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_4hfem")
tiles.4h.mal <- tileMethylCounts(myobj_4h.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_4hmal")

#unite calls
DMR.meth.4h <- unite(tiles.4h, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_4h")
DMR.meth.4h.fem <- unite(tiles.4h.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_4hfem")
DMR.meth.4h.mal <- unite(tiles.4h.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_4hmal")

#calculate differential methylation 
DMR.myDiff.4h <- calculateDiffMeth(DMR.meth.4h, mc.cores=1, covariates=covariates.4h, save.db = TRUE, suffix = "DMR_myDiff_4h")
DMR.myDiff.4h.fem <- calculateDiffMeth(DMR.meth.4h.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_4h")
DMR.myDiff.4h.mal <- calculateDiffMeth(DMR.meth.4h.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_4h")

#call significant methylation
DMR.diffMeth.4h <- getMethylDiff(DMR.myDiff.4h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_4h")
DMR.diffMeth.4h.fem <- getMethylDiff(DMR.myDiff.4h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_4hfem")
DMR.diffMeth.4h.mal <- getMethylDiff(DMR.myDiff.4h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_4hmal")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/05_shortterm_methylCalling_4hr-backup.RData")

q(save="yes")

