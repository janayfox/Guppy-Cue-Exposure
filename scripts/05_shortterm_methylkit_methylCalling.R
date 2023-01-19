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

setwd("/scratch/janayfox/guppyWGBS_shortterm/")
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

covariates.05h <- data.frame(tank=c("AC10","AC10","AC16","AC16","AC3","AC3",
                                    "C10","C10","C13","C13","C3","C3"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

covariates.1h <- data.frame(tank=c("AC14","AC14","AC2","AC2","AC9","AC9",
                                    "C14","C14","C3","C3","C9","C9"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

covariates.4h <- data.frame(tank=c("AC13","AC13","AC1","AC1","AC7","AC7",
                                    "C12","C12","C1","C1","C7","C7"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

covariates.24h <- data.frame(tank=c("AC15","AC15","AC4","AC4","AC8","AC8",
                                    "C15","C15","C4","C4","C8","C8"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

covariates.72h <- data.frame(tank=c("AC11","AC11","AC5","AC5","AC6","AC6",
                                    "C11","C11","C5","C5","C6","C6"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

## Find DMS ##
#remove sex chr (LG12) and unplaced scaffolds
myobj_all.subset <- selectByOverlap(myobj_all.filt, keep.chr.noXY)
myobj_all.fem.subset <- selectByOverlap(myobj_all.fem.filt, keep.chr.allchr)
myobj_all.mal.subset <- selectByOverlap(myobj_all.mal.filt, keep.chr.allchr)

myobj_05h.subset <- selectByOverlap(myobj_05h.filt, keep.chr.noXY)
myobj_05h.fem.subset <- selectByOverlap(myobj_05h.fem.filt, keep.chr.allchr)
myobj_05h.mal.subset <- selectByOverlap(myobj_05h.mal.filt, keep.chr.allchr)

myobj_1h.subset <- selectByOverlap(myobj_1h.filt, keep.chr.noXY)
myobj_1h.fem.subset <- selectByOverlap(myobj_1h.fem.filt, keep.chr.allchr)
myobj_1h.mal.subset <- selectByOverlap(myobj_1h.mal.filt, keep.chr.allchr)

myobj_4h.subset <- selectByOverlap(myobj_4h.filt, keep.chr.noXY)
myobj_4h.fem.subset <- selectByOverlap(myobj_4h.fem.filt, keep.chr.allchr)
myobj_4h.mal.subset <- selectByOverlap(myobj_4h.mal.filt, keep.chr.allchr)

myobj_24h.subset <- selectByOverlap(myobj_24h.filt, keep.chr.noXY)
myobj_24h.fem.subset <- selectByOverlap(myobj_24h.fem.filt, keep.chr.allchr)
myobj_24h.mal.subset <- selectByOverlap(myobj_24h.mal.filt, keep.chr.allchr)

myobj_72h.subset <- selectByOverlap(myobj_72h.filt, keep.chr.noXY)
myobj_72h.fem.subset <- selectByOverlap(myobj_72h.fem.filt, keep.chr.allchr)
myobj_72h.mal.subset <- selectByOverlap(myobj_72h.mal.filt, keep.chr.allchr)

#convert back to DB
myobj_all.subsetDB <- makeMethylDB(myobj_all.subset, "all_subset_DB")
myobj_all.fem.subsetDB <- makeMethylDB(myobj_all.fem.subset, "all_fem_subset_DB")
myobj_all.mal.subsetDB <- makeMethylDB(myobj_all.mal.subset, "all_mal_subset_DB")

myobj_05h.subsetDB <- makeMethylDB(myobj_05h.subset, "05h_subset_DB")
myobj_05h.fem.subsetDB <- makeMethylDB(myobj_05h.fem.subset, "05h_fem_subset_DB")
myobj_05h.mal.subsetDB <- makeMethylDB(myobj_05h.mal.subset, "05h_mal_subset_DB")

myobj_1h.subsetDB <- makeMethylDB(myobj_1h.subset, "1h_subset_DB")
myobj_1h.fem.subsetDB <- makeMethylDB(myobj_1h.fem.subset, "1h_fem_subset_DB")
myobj_1h.mal.subsetDB <- makeMethylDB(myobj_1h.mal.subset, "1h_mal_subset_DB")

myobj_4h.subsetDB <- makeMethylDB(myobj_4h.subset, "4h_subset_DB")
myobj_4h.fem.subsetDB <- makeMethylDB(myobj_4h.fem.subset, "4h_fem_subset_DB")
myobj_4h.mal.subsetDB <- makeMethylDB(myobj_4h.mal.subset, "4h_mal_subset_DB")

myobj_24h.subsetDB <- makeMethylDB(myobj_24h.subset, "24h_subset_DB")
myobj_24h.fem.subsetDB <- makeMethylDB(myobj_24h.fem.subset, "24h_fem_subset_DB")
myobj_24h.mal.subsetDB <- makeMethylDB(myobj_24h.mal.subset, "24h_mal_subset_DB")

myobj_72h.subsetDB <- makeMethylDB(myobj_72h.subset, "72h_subset_DB")
myobj_72h.fem.subsetDB <- makeMethylDB(myobj_72h.fem.subset, "72h_fem_subset_DB")
myobj_72h.mal.subsetDB <- makeMethylDB(myobj_72h.mal.subset, "72h_mal_subset_DB")

#filter out sites for 5x coverage 
myobj_all.5X=filterByCoverage(myobj_all.subsetDB,lo.count=5,lo.perc=NULL,
                                 hi.count=NULL, hi.perc=NULL, suffix = "all_5X")
myobj_all.fem.5X=filterByCoverage(myobj_all.fem.subsetDB,lo.count=5,lo.perc=NULL,
                              hi.count=NULL, hi.perc=NULL, suffix = "all_fem_5X")
myobj_all.mal.5X=filterByCoverage(myobj_all.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "all_mal_5X")

myobj_05h.5X=filterByCoverage(myobj_05h.subsetDB,lo.count=5,lo.perc=NULL,
                              hi.count=NULL, hi.perc=NULL, suffix = "05h_5X")
myobj_05h.fem.5X=filterByCoverage(myobj_05h.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "05h_fem_5X")
myobj_05h.mal.5X=filterByCoverage(myobj_05h.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "05h_mal_5X")

myobj_1h.5X=filterByCoverage(myobj_1h.subsetDB,lo.count=5,lo.perc=NULL,
                              hi.count=NULL, hi.perc=NULL, suffix = "1h_5X")
myobj_1h.fem.5X=filterByCoverage(myobj_1h.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "1h_fem_5X")
myobj_1h.mal.5X=filterByCoverage(myobj_1h.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "1h_mal_5X")

myobj_4h.5X=filterByCoverage(myobj_4h.subsetDB,lo.count=5,lo.perc=NULL,
                              hi.count=NULL, hi.perc=NULL, suffix = "4h_5X")
myobj_4h.fem.5X=filterByCoverage(myobj_4h.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "4h_fem_5X")
myobj_4h.mal.5X=filterByCoverage(myobj_4h.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "4h_mal_5X")

myobj_24h.5X=filterByCoverage(myobj_24h.subsetDB,lo.count=5,lo.perc=NULL,
                              hi.count=NULL, hi.perc=NULL, suffix = "24h_5X")
myobj_24h.fem.5X=filterByCoverage(myobj_24h.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "24h_fem_5X")
myobj_24h.mal.5X=filterByCoverage(myobj_24h.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "24h_mal_5X")

myobj_72h.5X=filterByCoverage(myobj_72h.subsetDB,lo.count=5,lo.perc=NULL,
                              hi.count=NULL, hi.perc=NULL, suffix = "72h_5X")
myobj_72h.fem.5X=filterByCoverage(myobj_72h.fem.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "72h_fem_5X")
myobj_72h.mal.5X=filterByCoverage(myobj_72h.mal.subsetDB,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL, hi.perc=NULL, suffix = "72h_mal_5X")

#normalize by median coverage
myobj_all.norm=normalizeCoverage(myobj_all.5X, method="all_median")
myobj_all.fem.norm=normalizeCoverage(myobj_all.fem.5X, method="all_fem_median")
myobj_all.mal.norm=normalizeCoverage(myobj_all.mal.5X, method="all_mal_median")

myobj_05h.norm=normalizeCoverage(myobj_05h.5X, method="05h_median")
myobj_05h.fem.norm=normalizeCoverage(myobj_05h.fem.5X, method="05h_fem_median")
myobj_05h.mal.norm=normalizeCoverage(myobj_05h.mal.5X, method="05h_mal_median")

myobj_1h.norm=normalizeCoverage(myobj_1h.5X, method="1h_median")
myobj_1h.fem.norm=normalizeCoverage(myobj_1h.fem.5X, method="1h_fem_median")
myobj_1h.mal.norm=normalizeCoverage(myobj_1h.mal.5X, method="1h_mal_median")

myobj_4h.norm=normalizeCoverage(myobj_4h.5X, method="4h_median")
myobj_4h.fem.norm=normalizeCoverage(myobj_4h.fem.5X, method="4h_fem_median")
myobj_4h.mal.norm=normalizeCoverage(myobj_4h.mal.5X, method="4h_mal_median")

myobj_24h.norm=normalizeCoverage(myobj_24h.5X, method="24h_median")
myobj_24h.fem.norm=normalizeCoverage(myobj_24h.fem.5X, method="24h_fem_median")
myobj_24h.mal.norm=normalizeCoverage(myobj_24h.mal.5X, method="24h_mal_median")

myobj_72h.norm=normalizeCoverage(myobj_72h.5X, method="72h_median")
myobj_72h.fem.norm=normalizeCoverage(myobj_72h.fem.5X, method="72h_fem_median")
myobj_72h.mal.norm=normalizeCoverage(myobj_72h.mal.5X, method="72h_mal_median")

#unite sites 
meth.all=unite(myobj_all.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_all")
meth.all.fem=unite(myobj_all.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_all_fem")
meth.all.mal=unite(myobj_all.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_all_mal")

meth.05h=unite(myobj_05h.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_05h")
meth.05h.fem=unite(myobj_05h.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_05h_fem")
meth.05h.mal=unite(myobj_05h.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_05h_mal")

meth.1h=unite(myobj_1h.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_1h")
meth.1h.fem=unite(myobj_1h.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_1h_fem")
meth.1h.mal=unite(myobj_1h.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_1h_mal")

meth.4h=unite(myobj_1h.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_4h")
meth.4h.fem=unite(myobj_1h.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_4h_fem")
meth.4h.mal=unite(myobj_1h.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_4h_mal")

meth.24h=unite(myobj_1h.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_24h")
meth.24h.fem=unite(myobj_1h.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_24h_fem")
meth.24h.mal=unite(myobj_1h.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_24h_mal")

meth.72h=unite(myobj_1h.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_72h")
meth.72h.fem=unite(myobj_1h.fem.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_72h_fem")
meth.72h.mal=unite(myobj_1h.mal.norm, destrand=TRUE, save.db = TRUE, suffix = "DMS_unite_72h_mal")

#calculate differential methylation
DMS.myDiff.all <- calculateDiffMeth(meth.all, covariates=covariates.all, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_all")
DMS.myDiff.all.fem <- calculateDiffMeth(meth.all.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_all_fem")
DMS.myDiff.all.mal <- calculateDiffMeth(meth.all.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_all_mal")

DMS.myDiff.05h <- calculateDiffMeth(meth.05h, covariates=covariates.05h, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_05h")
DMS.myDiff.05h.fem <- calculateDiffMeth(meth.05h.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_05h_fem")
DMS.myDiff.05h.mal <- calculateDiffMeth(meth.05h.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_05h_mal")

DMS.myDiff.1h <- calculateDiffMeth(meth.1h, covariates=covariates.1h, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_1h")
DMS.myDiff.1h.fem <- calculateDiffMeth(meth.1h.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_1h_fem")
DMS.myDiff.1h.mal <- calculateDiffMeth(meth.1h.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_1h_mal")

DMS.myDiff.4h <- calculateDiffMeth(meth.4h, covariates=covariates.4h, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_4h")
DMS.myDiff.4h.fem <- calculateDiffMeth(meth.4h.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_4h_fem")
DMS.myDiff.4h.mal <- calculateDiffMeth(meth.4h.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_4h_mal")

DMS.myDiff.24h <- calculateDiffMeth(meth.24h, covariates=covariates.24h, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_24h")
DMS.myDiff.24h.fem <- calculateDiffMeth(meth.24h.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_24h_fem")
DMS.myDiff.24h.mal <- calculateDiffMeth(meth.24h.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_24h_mal")

DMS.myDiff.72h <- calculateDiffMeth(meth.72h, covariates=covariates.72h, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_72h")
DMS.myDiff.72h.fem <- calculateDiffMeth(meth.72h.fem, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_72h_fem")
DMS.myDiff.72h.mal <- calculateDiffMeth(meth.72h.mal, mc.cores=1, save.db = TRUE, suffix = "DMS_myDiff_72h_mal")

#call significant methylation
DMS.diffMeth.all <- getMethylDiff(DMS.myDiff.all, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_all")
DMS.diffMeth.all.fem <- getMethylDiff(DMS.myDiff.all.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_all.fem")
DMS.diffMeth.all.mal <- getMethylDiff(DMS.myDiff.all.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_all.mal")

DMS.diffMeth.05h <- getMethylDiff(DMS.myDiff.05h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_05h")
DMS.diffMeth.05h.fem <- getMethylDiff(DMS.myDiff.05h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_05h.fem")
DMS.diffMeth.05h.mal <- getMethylDiff(DMS.myDiff.05h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_05h.mal")

DMS.diffMeth.1h <- getMethylDiff(DMS.myDiff.1h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_1h")
DMS.diffMeth.1h.fem <- getMethylDiff(DMS.myDiff.1h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_1h.fem")
DMS.diffMeth.1h.mal <- getMethylDiff(DMS.myDiff.1h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_1h.mal")

DMS.diffMeth.4h <- getMethylDiff(DMS.myDiff.4h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_4h")
DMS.diffMeth.4h.fem <- getMethylDiff(DMS.myDiff.4h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_4h.fem")
DMS.diffMeth.4h.mal <- getMethylDiff(DMS.myDiff.4h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_4h.mal")

DMS.diffMeth.24h <- getMethylDiff(DMS.myDiff.24h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_24h")
DMS.diffMeth.24h.fem <- getMethylDiff(DMS.myDiff.24h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_24h.fem")
DMS.diffMeth.24h.mal <- getMethylDiff(DMS.myDiff.24h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_24h.mal")

DMS.diffMeth.72h <- getMethylDiff(DMS.myDiff.72h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_72h")
DMS.diffMeth.72h.fem <- getMethylDiff(DMS.myDiff.72h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_72h.fem")
DMS.diffMeth.72h.mal <- getMethylDiff(DMS.myDiff.72h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_72h.mal")

## Find DMRs ##
#tile into 100 bp windows with min coverage 10X and 5X
tiles.all <- tileMethylCounts(myobj_all.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_all")
tiles.all.fem <- tileMethylCounts(myobj_all.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_allfem")
tiles.all.mal <- tileMethylCounts(myobj_all.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_allmal")

tiles.05h <- tileMethylCounts(myobj_05h.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_05h")
tiles.05h.fem <- tileMethylCounts(myobj_05h.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_05hfem")
tiles.05h.mal <- tileMethylCounts(myobj_05h.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_05hmal")

tiles.1h <- tileMethylCounts(myobj_1h.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_1h")
tiles.1h.fem <- tileMethylCounts(myobj_1h.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_1hfem")
tiles.1h.mal <- tileMethylCounts(myobj_1h.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_1hmal")

tiles.4h <- tileMethylCounts(myobj_4h.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_4h")
tiles.4h.fem <- tileMethylCounts(myobj_4h.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_4hfem")
tiles.4h.mal <- tileMethylCounts(myobj_4h.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_4hmal")

tiles.24h <- tileMethylCounts(myobj_24h.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_24h")
tiles.24h.fem <- tileMethylCounts(myobj_24h.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_24hfem")
tiles.24h.mal <- tileMethylCounts(myobj_24h.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_24hmal")

tiles.72h <- tileMethylCounts(myobj_72h.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_72h")
tiles.72h.fem <- tileMethylCounts(myobj_72h.fem.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_72hfem")
tiles.72h.mal <- tileMethylCounts(myobj_72h.mal.subsetDB, win.size = 100, step.size = 100, cov.bases = 5, save.db = TRUE, suffix = "tiles_72hmal")

#unite calls
DMR.meth.all <- unite(tiles.all, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_all")
DMR.meth.all.fem <- unite(tiles.all.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_allfem")
DMR.meth.all.mal <- unite(tiles.all.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_allmal")

DMR.meth.05h <- unite(tiles.05h, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_05h")
DMR.meth.05h.fem <- unite(tiles.05h.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_05hfem")
DMR.meth.05h.mal <- unite(tiles.05h.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_05hmal")

DMR.meth.1h <- unite(tiles.1h, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_1h")
DMR.meth.1h.fem <- unite(tiles.1h.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_1hfem")
DMR.meth.1h.mal <- unite(tiles.1h.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_1hmal")

DMR.meth.4h <- unite(tiles.4h, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_4h")
DMR.meth.4h.fem <- unite(tiles.4h.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_4hfem")
DMR.meth.4h.mal <- unite(tiles.4h.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_4hmal")

DMR.meth.24h <- unite(tiles.24h, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_24h")
DMR.meth.24h.fem <- unite(tiles.24h.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_24hfem")
DMR.meth.24h.mal <- unite(tiles.24h.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_24hmal")

DMR.meth.72h <- unite(tiles.72h, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_72h")
DMR.meth.72h.fem <- unite(tiles.72h.fem, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_72hfem")
DMR.meth.72h.mal <- unite(tiles.72h.mal, destrand=TRUE, save.db = TRUE, suffix = "DMR_unite_72hmal")

#calculate differential methylation 
DMR.myDiff.all <- calculateDiffMeth(DMR.meth.all, mc.cores=1, covariates=covariates.all, save.db = TRUE, suffix = "DMR_myDiff_all")
DMR.myDiff.all.fem <- calculateDiffMeth(DMR.meth.all.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_all")
DMR.myDiff.all.mal <- calculateDiffMeth(DMR.meth.all.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_all")

DMR.myDiff.05h <- calculateDiffMeth(DMR.meth.05h, mc.cores=1, covariates=covariates.05h, save.db = TRUE, suffix = "DMR_myDiff_05h")
DMR.myDiff.05h.fem <- calculateDiffMeth(DMR.meth.05h.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_05h")
DMR.myDiff.05h.mal <- calculateDiffMeth(DMR.meth.05h.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_05h")

DMR.myDiff.1h <- calculateDiffMeth(DMR.meth.1h, mc.cores=1, covariates=covariates.1h, save.db = TRUE, suffix = "DMR_myDiff_1h")
DMR.myDiff.1h.fem <- calculateDiffMeth(DMR.meth.1h.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_1h")
DMR.myDiff.1h.mal <- calculateDiffMeth(DMR.meth.1h.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_1h")

DMR.myDiff.4h <- calculateDiffMeth(DMR.meth.4h, mc.cores=1, covariates=covariates.4h, save.db = TRUE, suffix = "DMR_myDiff_4h")
DMR.myDiff.4h.fem <- calculateDiffMeth(DMR.meth.4h.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_4h")
DMR.myDiff.4h.mal <- calculateDiffMeth(DMR.meth.4h.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_4h")

DMR.myDiff.24h <- calculateDiffMeth(DMR.meth.24h, mc.cores=1, covariates=covariates.24h, save.db = TRUE, suffix = "DMR_myDiff_24h")
DMR.myDiff.24h.fem <- calculateDiffMeth(DMR.meth.24h.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_24h")
DMR.myDiff.24h.mal <- calculateDiffMeth(DMR.meth.24h.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_24h")

DMR.myDiff.72h <- calculateDiffMeth(DMR.meth.72h, mc.cores=1, covariates=covariates.72h, save.db = TRUE, suffix = "DMR_myDiff_72h")
DMR.myDiff.72h.fem <- calculateDiffMeth(DMR.meth.72h.fem, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_72h")
DMR.myDiff.72h.mal <- calculateDiffMeth(DMR.meth.72h.mal, mc.cores=1, save.db = TRUE, suffix = "DMR_myDiff_72h")

#call significant methylation
DMR.diffMeth.all <- getMethylDiff(DMR.myDiff.all, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_all")
DMR.diffMeth.all.fem <- getMethylDiff(DMR.myDiff.all.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_allfem")
DMR.diffMeth.all.mal <- getMethylDiff(DMR.myDiff.all.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_allmal")

DMR.diffMeth.05h <- getMethylDiff(DMR.myDiff.05h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05h")
DMR.diffMeth.05h.fem <- getMethylDiff(DMR.myDiff.05h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05hfem")
DMR.diffMeth.05h.mal <- getMethylDiff(DMR.myDiff.05h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05hmal")

DMR.diffMeth.1h <- getMethylDiff(DMR.myDiff.1h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_1h")
DMR.diffMeth.1h.fem <- getMethylDiff(DMR.myDiff.1h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_1hfem")
DMR.diffMeth.1h.mal <- getMethylDiff(DMR.myDiff.1h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_1hmal")

DMR.diffMeth.4h <- getMethylDiff(DMR.myDiff.4h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_4h")
DMR.diffMeth.4h.fem <- getMethylDiff(DMR.myDiff.4h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_4hfem")
DMR.diffMeth.4h.mal <- getMethylDiff(DMR.myDiff.4h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_4hmal")

DMR.diffMeth.24h <- getMethylDiff(DMR.myDiff.24h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_24h")
DMR.diffMeth.24h.fem <- getMethylDiff(DMR.myDiff.24h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_24hfem")
DMR.diffMeth.24h.mal <- getMethylDiff(DMR.myDiff.24h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_24hmal")

DMR.diffMeth.72h <- getMethylDiff(DMR.myDiff.72h, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_72h")
DMR.diffMeth.72h.fem <- getMethylDiff(DMR.myDiff.72h.fem, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_72hfem")
DMR.diffMeth.72h.mal <- getMethylDiff(DMR.myDiff.72h.mal, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_72hmal")

## Save workspace image for later loading ##
save.image(file = ".RData")
save.image(file = "./backupRData/05_shortterm_methylCalling-backup.RData")

q(save="yes")

