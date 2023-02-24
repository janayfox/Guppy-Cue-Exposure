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

setwd("/scratch/janayfox/guppyWGBS_shortterm/72hr/")
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
covariates.72h <- data.frame(tank=c("AC11","AC11","AC5","AC5","AC6","AC6",
                                    "C11","C11","C5","C5","C6","C6"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

## Find DMS ##
#remove sex chr (LG12) and unplaced scaffolds
myobj_72h.subset <- selectByOverlap(myobj_72h.filt, keep.chr.noXY)
myobj_72h.fem.subset <- selectByOverlap(myobj_72h.fem.filt, keep.chr.allchr)
myobj_72h.mal.subset <- selectByOverlap(myobj_72h.mal.filt, keep.chr.allchr)

#convert back to DB
myobj_72h.subsetDB <- makeMethylDB(myobj_72h.subset, "72h_subset_DB")
myobj_72h.fem.subsetDB <- makeMethylDB(myobj_72h.fem.subset, "72h_fem_subset_DB")
myobj_72h.mal.subsetDB <- makeMethylDB(myobj_72h.mal.subset, "72h_mal_subset_DB")

#filter out sites for 5x coverage 
myobj_72h.5X=filterByCoverage(myobj_72h.subsetDB,lo.count=5,lo.perc=NULL,
                              hi.count=NULL, hi.perc=NULL, suffix = "72h_5X")
myobj_72h.fem.5X=filterByCoverage(myobj_72h.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "72h_fem_5X")
myobj_72h.mal.5X=filterByCoverage(myobj_72h.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "72h_mal_5X")

#normalize by median coverage
myobj_72h.norm=normalizeCoverage(myobj_72h.5X, method="median")
myobj_72h.fem.norm=normalizeCoverage(myobj_72h.fem.5X, method="median")
myobj_72h.mal.norm=normalizeCoverage(myobj_72h.mal.5X, method="median")

#unite sites 
meth.72h=unite(myobj_72h.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_72h")
meth.72h.fem=unite(myobj_72h.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_72h_fem")
meth.72h.mal=unite(myobj_72h.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_72h_mal")

#calculate differential methylation
DMS.myDiff.72h <- calculateDiffMeth(meth.72h, covariates=covariates.72h, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_72h")
DMS.myDiff.72h.fem <- calculateDiffMeth(meth.72h.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_72h_fem")
DMS.myDiff.72h.mal <- calculateDiffMeth(meth.72h.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_72h_mal")

#call significant methylation
DMS.diffMeth.72h <- getMethylDiff(DMS.myDiff.72h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_72h")
DMS.diffMeth.72h.fem <- getMethylDiff(DMS.myDiff.72h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_72h.fem")
DMS.diffMeth.72h.mal <- getMethylDiff(DMS.myDiff.72h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_72h.mal")

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X and 5X
tiles.72h <- tileMethylCounts(myobj_72h.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_72h")
tiles.72h.fem <- tileMethylCounts(myobj_72h.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_72hfem")
tiles.72h.mal <- tileMethylCounts(myobj_72h.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_72hmal")

#unite calls
DMR.meth.72h <- unite(tiles.72h, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_72h")
DMR.meth.72h.fem <- unite(tiles.72h.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_72hfem")
DMR.meth.72h.mal <- unite(tiles.72h.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_72hmal")

#calculate differential methylation 
DMR.myDiff.72h <- calculateDiffMeth(DMR.meth.72h, mc.cores=1, covariates=covariates.72h, save.db = TRUE, suffix = "DMR_myDiff_72h")
DMR.myDiff.72h.fem <- calculateDiffMeth(DMR.meth.72h.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_72h")
DMR.myDiff.72h.mal <- calculateDiffMeth(DMR.meth.72h.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_72h")

#call significant methylation
DMR.diffMeth.72h <- getMethylDiff(DMR.myDiff.72h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_72h")
DMR.diffMeth.72h.fem <- getMethylDiff(DMR.myDiff.72h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_72hfem")
DMR.diffMeth.72h.mal <- getMethylDiff(DMR.myDiff.72h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_72hmal")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/05_shortterm_methylCalling_72hr-backup.RData")

q(save="yes")

