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
library("genomation", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/methylKit/st/05h")

## Prepare tabix files
#create lists of file locations
file.list.05h = list("../../../mergedCov/st/ST2AC10F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2AC10M.CpG_merged.cov",
                     "../../../mergedCov/st/ST2AC16F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2AC16M.CpG_merged.cov",
                     "../../../mergedCov/st/ST2AC3F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2AC3M.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C10F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C10M.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C13F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C13M.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C2F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C2M.CpG_merged.cov")

file.list.05h.fem = list("../../../mergedCov/st/ST2AC10F.CpG_merged.cov",
                         "../../../mergedCov/st/ST2AC16F.CpG_merged.cov",
                         "../../../mergedCov/st/ST2AC3F.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C10F.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C13F.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C2F.CpG_merged.cov")

file.list.05h.mal = list("../../../mergedCov/st/ST2AC10M.CpG_merged.cov",
                         "../../../mergedCov/st/ST2AC16M.CpG_merged.cov",
                         "../../../mergedCov/st/ST2AC3M.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C10M.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C13M.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C2M.CpG_merged.cov")

#create tabix file
myobj.05h=methRead(file.list.05h,
                   sample.id=list("ST2AC10F","ST2AC10M","ST2AC16F","ST2AC16M","ST2AC3F","ST2AC3M",
                                  "ST2C10F","ST2C10M","ST2C13F","ST2C13M","ST2C3F","ST2C3M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCoverage",
                   treatment=c(1,1,1,1,1,1,
                               0,0,0,0,0,0),
                   context="CpG",
                   mincov = 1,
                   dbtype = "tabix",
                   dbdir = "shortterm_05h_DB_merged"
)

myobj.05h.fem=methRead(file.list.05h.fem,
                       sample.id=list("ST2AC10F","ST2AC16F","ST2AC3F",
                                      "ST2C10F","ST2C13F","ST2C3F"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,
                                   0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_05hF_DB_merged"
)

myobj.05h.mal=methRead(file.list.05h.mal,
                       sample.id=list("ST2AC10M","ST2AC16M","ST2AC3M",
                                      "ST2C10M","ST2C13M","ST2C3M"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,
                                   0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_05hM_DB_merged"
)

#get coverage stats 
getCoverageStats(myobj.05h[[2]], both.strands = FALSE)
getCoverageStats(myobj.05h.fem[[2]], both.strands = FALSE)
getCoverageStats(myobj.05h.mal[[2]], both.strands = FALSE)

#filter out sites in the 99.9th percentile of coverage (PCR bias) 
myobj.05h.5X=filterByCoverage(myobj.05h,lo.count=5,lo.perc=NULL,
                               hi.count=NULL, hi.perc=99.9, suffix = "05h_5X_merged")
myobj.05h.fem.5X=filterByCoverage(myobj.05h.fem,lo.count=5,lo.perc=NULL,
                                   hi.count=NULL, hi.perc=99.9, suffix = "05h_5X_fem_merged")
myobj.05h.mal.5X=filterByCoverage(myobj.05h.mal,lo.count=5,lo.perc=NULL,
                                   hi.count=NULL, hi.perc=99.9, suffix = "05h_5X_mal_merged")

#normalize by median coverage
norm.myobj.05h.5X=normalizeCoverage(myobj.05h.5X, method="median")
norm.myobj.05h.fem.5X=normalizeCoverage(myobj.05h.fem.5X, method="median")
norm.myobj.05h.mal.5X=normalizeCoverage(myobj.05h.mal.5X, method="median")

## Remove sex chromosomes and unplacex scaffolds ##
#prepare GRanges object for chromosomes to keep 
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

#remove sex chromosomes and unplaced scaffolds
myobj.05h.subset <- selectByOverlap(norm.myobj.05h.5X, keep.chr.noXY)
myobj.05h.fem.subset <- selectByOverlap(norm.myobj.05h.fem.5X, keep.chr.allchr)
myobj.05h.mal.subset <- selectByOverlap(norm.myobj.05h.mal.5X, keep.chr.allchr)

## Find DMRs ##
#enter covariates 
covariates.05h <- data.frame(tank=c("AC10","AC10","AC16","AC16","AC3","AC3",
                                    "C10","C10","C13","C13","C3","C3"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

#tile into 100 bp windows with min coverage 10X and 5X
tiles.05h.10X <- tileMethylCounts(myobj.05h.subset, win.size = 100, step.size = 100, cov.bases = 10)
tiles.05h.fem.10X <- tileMethylCounts(myobj.05h.fem.subset, win.size = 100, step.size = 100, cov.bases = 10)
tiles.05h.mal.10X <- tileMethylCounts(myobj.05h.mal.subset, win.size = 100, step.size = 100, cov.bases = 10)

tiles.05h.5X <- tileMethylCounts(myobj.05h.subset, win.size = 100, step.size = 100, cov.bases = 5)
tiles.05h.fem.5X <- tileMethylCounts(myobj.05h.fem.subset, win.size = 100, step.size = 100, cov.bases = 5)
tiles.05h.mal.5X <- tileMethylCounts(myobj.05h.mal.subset, win.size = 100, step.size = 100, cov.bases = 5)

#check number of tiles 
tiles.05h.10X
tiles.05h.fem.10X
tiles.05h.mal.10X

tiles.05h.5X
tiles.05h.fem.5X
tiles.05h.mal.5X

#unite calls
DMR.meth.05h.10X <- unite(tiles.05h.10X, min.per.group = 6L, destrand=FALSE, save.db = TRUE, suffix = "DMR_unite_05h_10X")
DMR.meth.05h.fem.10X <- unite(tiles.05h.fem.10X, min.per.group = myobj.3XL, destrand=FALSE, save.db = TRUE, suffix = "DMR_unite_05h_fem_10X")
DMR.meth.05h.mal.10X <- unite(tiles.05h.mal.10X, min.per.group = 3L, destrand=FALSE, save.db = TRUE, suffix = "DMR_unite_05h_mal_10X")

DMR.meth.05h.5X <- unite(tiles.05h.5X, destrand=FALSE, min.per.group = 6L, save.db = TRUE, suffix = "DMR_unite_05h_5X")
DMR.meth.05h.fem.5X <- unite(tiles.05h.fem.5X, destrand=FALSE, min.per.group = 3L, save.db = TRUE, suffix = "DMR_unite_05h_fem_5X")
DMR.meth.05h.mal.5X <- unite(tiles.05h.mal.5X, destrand=FALSE, min.per.group = 3L, save.db = TRUE, suffix = "DMR_unite_05h_mal_5X")

#check number of regions retained 
DMR.meth.05h.10X
DMR.meth.05h.fem.10X
DMR.meth.05h.mal.10X

DMR.meth.05h.5X
DMR.meth.05h.fem.5X
DMR.meth.05h.mal.5X

#calculate differential methylation 
DMR.myDiff.05h.10X <- calculateDiffMeth(DMR.meth.05h.10X, mc.cores=2, test="Chisq", covariates=covariates.05h, save.db = TRUE, suffix = "DMR_myDiff_05h_10X")
DMR.myDiff.05h.fem.10X <- calculateDiffMeth(DMR.meth.05h.fem.10X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMR_myDiff_05h_fem_10X")
DMR.myDiff.05h.mal.10X <- calculateDiffMeth(DMR.meth.05h.mal.10X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMR_myDiff_05h_mal_10X")

DMR.myDiff.05h.5X <- calculateDiffMeth(DMR.meth.05h.5X, mc.cores=2, test="Chisq", covariates=covariates.05h, save.db = TRUE, suffix = "DMR_myDiff_05h_5X")
DMR.myDiff.05h.fem.5X <- calculateDiffMeth(DMR.meth.05h.fem.5X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMR_myDiff_05h_fem_5X")
DMR.myDiff.05h.mal.5X <- calculateDiffMeth(DMR.meth.05h.mal.5X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMR_myDiff_05h_mal_5X")

#call significant methylation
DMR.diffMeth.05h.10X <- getMethylDiff(DMR.myDiff.05h.10X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05h_10X")
DMR.diffMeth.05h.fem.10X <- getMethylDiff(DMR.myDiff.05h.fem.10X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05h_fem_10X")
DMR.diffMeth.05h.mal.10X <- getMethylDiff(DMR.myDiff.05h.mal.10X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05h_mal_10X")

DMR.diffMeth.05h.5X <- getMethylDiff(DMR.myDiff.05h.5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05h_5X")
DMR.diffMeth.05h.fem.5X <- getMethylDiff(DMR.myDiff.05h.fem.5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05h_fem_5X")
DMR.diffMeth.05h.mal.5X <- getMethylDiff(DMR.myDiff.05h.mal.5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMR_diffMeth_05h_mal_5X")

#check number of DMRs 
DMR.diffMeth.05h.10X
DMR.diffMeth.05h.fem.10X
DMR.diffMeth.05h.mal.10X

DMR.diffMeth.05h.5X
DMR.diffMeth.05h.fem.5X
DMR.diffMeth.05h.mal.5X

#get meth per chromosome 
DMR.diffMeth.05h.10X.chr <- diffMethPerChr(DMR.diffMeth.05h.10X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "chrDMR_05h_10X")
DMR.diffMeth.05h.10X.chr
DMR.diffMeth.05h.fem.10X.chr <- diffMethPerChr(DMR.diffMeth.05h.fem.10X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "chrDMR_05h_fem_10X")
DMR.diffMeth.05h.fem.10X.chr
DMR.diffMeth.05h.mal.10X.chr <- diffMethPerChr(DMR.diffMeth.05h.mal.10X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "chrDMR_05h_mal_10X")
DMR.diffMeth.05h.mal.10X.chr

DMR.diffMeth.05h.5X.chr <- diffMethPerChr(DMR.diffMeth.05h.5X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "chrDMR_05h_5X")
DMR.diffMeth.05h.5X.chr
DMR.diffMeth.05h.fem.5X.chr <- diffMethPerChr(DMR.diffMeth.05h.fem.5X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "chrDMR_05h_fem_5X")
DMR.diffMeth.05h.fem.5X.chr
DMR.diffMeth.05h.mal.5X.chr <- diffMethPerChr(DMR.diffMeth.05h.mal.5X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "chrDMR_05h_mal_5X")
DMR.diffMeth.05h.mal.5X.chr

## Save R objects ##
saveRDS(DMR.meth.05h.10X, file = "./DMRmeth_05h_10X.RDS")
saveRDS(DMR.meth.05h.fem.10X, file = "./DMRmeth_05h_fem_10X.RDS")
saveRDS(DMR.meth.05h.mal.10X, file = "./DMRmeth_05h_mal_10X.RDS")
saveRDS(DMR.meth.05h.5X, file = "./DMRmeth_05h_5X.RDS")
saveRDS(DMR.meth.05h.fem.5X, file = "./DMRmeth_05h_fem_5X.RDS")
saveRDS(DMR.meth.05h.mal.5X, file = "./DMRmeth_05h_mal_5X.RDS")

saveRDS(getData(DMR.meth.05h.10X), file = "./DMRmeth_05h_10X_data.RDS")
saveRDS(getData(DMR.meth.05h.fem.10X), file = "./DMRmeth_05h_fem_10X_data.RDS")
saveRDS(getData(DMR.meth.05h.mal.10X), file = "./DMRmeth_05h_mal_10X_data.RDS")
saveRDS(getData(DMR.meth.05h.5X), file = "./DMRmeth_05h_5X_data.RDS")
saveRDS(getData(DMR.meth.05h.fem.5X), file = "./DMRmeth_05h_fem_5X_data.RDS")
saveRDS(getData(DMR.meth.05h.mal.5X), file = "./DMRmeth_05h_mal_5X_data.RDS")

saveRDS(DMR.myDiff.05h.10X, file = "./DMRmyDiff_05h_10X.RDS")
saveRDS(DMR.myDiff.05h.fem.10X, file = "./DMRmyDiff_05h_fem_10X.RDS")
saveRDS(DMR.myDiff.05h.mal.10X, file = "./DMRmyDiff_05h_mal_10X.RDS")
saveRDS(DMR.myDiff.05h.5X, file = "./DMRmyDiff_05h_5X.RDS")
saveRDS(DMR.myDiff.05h.fem.5X, file = "./DMRmyDiff_05h_fem_5X.RDS")
saveRDS(DMR.myDiff.05h.mal.5X, file = "./DMRmyDiff_05h_mal_5X.RDS")

saveRDS(getData(DMR.myDiff.05h.10X), file = "./DMRmyDiff_05h_10X_data.RDS")
saveRDS(getData(DMR.myDiff.05h.fem.10X), file = "./DMRmyDiff_05h_fem_10X_data.RDS")
saveRDS(getData(DMR.myDiff.05h.mal.10X), file = "./DMRmyDiff_05h_mal_10X_data.RDS")
saveRDS(getData(DMR.myDiff.05h.5X), file = "./DMRmyDiff_05h_5X_data.RDS")
saveRDS(getData(DMR.myDiff.05h.fem.5X), file = "./DMRmyDiff_05h_fem_5X_data.RDS")
saveRDS(getData(DMR.myDiff.05h.mal.5X), file = "./DMRmyDiff_05h_mal_5X_data.RDS")

saveRDS(DMR.diffMeth.05h.10X, file = "./DMRdiffMeth_05h_10X.RDS")
saveRDS(DMR.diffMeth.05h.fem.10X, file = "./DMRdiffMeth_05h_fem_10X.RDS")
saveRDS(DMR.diffMeth.05h.mal.10X, file = "./DMRdiffMeth_05h_mal_10X.RDS")
saveRDS(DMR.diffMeth.05h.5X, file = "./DMRdiffMeth_05h_5X.RDS")
saveRDS(DMR.diffMeth.05h.fem.5X, file = "./DMRdiffMeth_05h_fem_5X.RDS")
saveRDS(DMR.diffMeth.05h.mal.5X, file = "./DMRdiffMeth_05h_mal_5X.RDS")

saveRDS(getData(DMR.diffMeth.05h.10X), file = "./DMRdiffMeth_05h_10X_data.RDS")
saveRDS(getData(DMR.diffMeth.05h.fem.10X), file = "./DMRdiffMeth_05h_fem_10X_data.RDS")
saveRDS(getData(DMR.diffMeth.05h.mal.10X), file = "./DMRdiffMeth_05h_mal_10X_data.RDS")
saveRDS(getData(DMR.diffMeth.05h.5X), file = "./DMRdiffMeth_05h_5X_data.RDS")
saveRDS(getData(DMR.diffMeth.05h.fem.5X), file = "./DMRdiffMeth_05h_fem_5X_data.RDS")
saveRDS(getData(DMR.diffMeth.05h.mal.5X), file = "./DMRdiffMeth_05h_mal_5X_data.RDS")

saveRDS(DMR.diffMeth.05h.10X.chr, file = "./DMRdiffMethChr_05h_10X.RDS")
saveRDS(DMR.diffMeth.05h.fem.10X.chr, file = "./DMRdiffMethChr_05h_fem_10X.RDS")
saveRDS(DMR.diffMeth.05h.mal.10X.chr, file = "./DMRdiffMethChr_05h_mal_10X.RDS")
saveRDS(DMR.diffMeth.05h.5X.chr, file = "./DMRdiffMethChr_05h_5X.RDS")
saveRDS(DMR.diffMeth.05h.fem.5X.chr, file = "./DMRdiffMethChr_05h_fem_5X.RDS")
saveRDS(DMR.diffMeth.05h.mal.5X.chr, file = "./DMRdiffMethChr_05h_mal_5X.RDS")
