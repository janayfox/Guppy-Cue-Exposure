#####################################################################################################################
### Goal: Read alignment files into methylKit and filter to create tabix files of filtered cytosine methylation
### Author: Janay Fox
### R script
#####################################################################################################################

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

setwd("/scratch/janayfox/guppyWGBS/methylKit/st/72h")

## Prepare tabix files
#create lists of file locations
file.list.72h = list("../../../mergedCov/st/ST2AC11F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2AC11M.CpG_merged.cov",
                     "../../../mergedCov/st/ST2AC5F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2AC5M.CpG_merged.cov",
                     "../../../mergedCov/st/ST2AC6F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2AC6M.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C11F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C11M.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C5F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C5M.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C6F.CpG_merged.cov",
                     "../../../mergedCov/st/ST2C6M.CpG_merged.cov")

file.list.72h.fem = list("../../../mergedCov/st/ST2AC11F.CpG_merged.cov",
                         "../../../mergedCov/st/ST2AC5F.CpG_merged.cov",
                         "../../../mergedCov/st/ST2AC6F.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C11F.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C5F.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C6F.CpG_merged.cov")

file.list.72h.mal = list("../../../mergedCov/st/ST2AC11M.CpG_merged.cov",
                         "../../../mergedCov/st/ST2AC5M.CpG_merged.cov",
                         "../../../mergedCov/st/ST2AC6M.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C11M.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C5M.CpG_merged.cov",
                         "../../../mergedCov/st/ST2C6M.CpG_merged.cov")

#create tabix file
myobj.72h=methRead(file.list.72h,
                   sample.id=list("ST2AC11F","ST2AC11M","ST2AC5F","ST2AC5M","ST2AC6F","ST2AC6M",
                                  "ST2C11F","ST2C11M","ST2C5F","ST2C5M","ST2C6F","ST2C6M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCoverage",
                   treatment=c(1,1,1,1,1,1,
                               0,0,0,0,0,0),
                   context="CpG",
                   mincov = 1,
                   dbtype = "tabix",
                   dbdir = "shortterm_72h_DB_merged"
)

myobj.72h.fem=methRead(file.list.72h.fem,
                       sample.id=list("ST2AC11F","ST2AC5F","ST2AC6F",
                                      "ST2C11F","ST2C5F","ST2C6F"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,
                                   0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_72hF_DB_merged"
)

myobj.72h.mal=methRead(file.list.72h.mal,
                       sample.id=list("ST2AC11M","ST2AC5M","ST2AC6M",
                                      "ST2C11M","ST2C5M","ST2C6M"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,
                                   0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_72hM_DB_merged"
)

#get coverage stats 
getCoverageStats(myobj.72h[[2]], both.strands = FALSE)
getCoverageStats(myobj.72h.fem[[2]], both.strands = FALSE)
getCoverageStats(myobj.72h.mal[[2]], both.strands = FALSE)

#filter out sites in the 99.9th percentile of coverage (PCR bias) 
myobj.72h.5X=filterByCoverage(myobj.72h,lo.count=5,lo.perc=NULL,
                               hi.count=NULL, hi.perc=99.9, suffix = "72h_5X_merged")
myobj.72h.fem.5X=filterByCoverage(myobj.72h.fem,lo.count=5,lo.perc=NULL,
                                   hi.count=NULL, hi.perc=99.9, suffix = "72h_5X_fem_merged")
myobj.72h.mal.5X=filterByCoverage(myobj.72h.mal,lo.count=5,lo.perc=NULL,
                                   hi.count=NULL, hi.perc=99.9, suffix = "72h_5X_mal_merged")

#normalize by median coverage
norm.myobj.72h.5X=normalizeCoverage(myobj.72h.5X, method="median")
norm.myobj.72h.fem.5X=normalizeCoverage(myobj.72h.fem.5X, method="median")
norm.myobj.72h.mal.5X=normalizeCoverage(myobj.72h.mal.5X, method="median")

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

#remove sex chr (LG12) and unplaced scaffolds
myobj.72h.subset <- selectByOverlap(norm.myobj.72h.5X, keep.chr.noXY)
myobj.72h.fem.subset <- selectByOverlap(norm.myobj.72h.fem.5X, keep.chr.allchr)
myobj.72h.mal.subset <- selectByOverlap(norm.myobj.72h.mal.5X, keep.chr.allchr)

##Find DMS##
#unite sites 
DMS.meth.72h.5X=unite(myobj.72h.subset, min.per.group = 6L, destrand=FALSE, save.db = TRUE, suffix = "DMS_unite_72h")
DMS.meth.72h.fem.5X=unite(myobj.72h.fem.subset, min.per.group = 3L, destrand=FALSE, save.db = TRUE, suffix = "DMS_unite_72h_fem")
DMS.meth.72h.mal.5X=unite(myobj.72h.mal.subset, min.per.group = 3L, destrand=FALSE, save.db = TRUE, suffix = "DMS_unite_72h_mal")

#filter out low variation sites 
pm.72h.5x <- percMethylation(DMS.meth.72h.5X) #get percent methylation matrix
sds.72h.5x <- matrixStats::rowSds(pm.72h.5x) #calculate standard deviation of CpGs 
DMS.meth.72h.5X <- DMS.meth.72h.5X[sds.72h.5x > 2]

pm.72h.fem.5x <- percMethylation(DMS.meth.72h.fem.5X) #get percent methylation matrix
sds.72h.fem.5x <- matrixStats::rowSds(pm.72h.fem.5x) #calculate standard deviation of CpGs 
DMS.meth.72h.fem.5X <- DMS.meth.72h.fem.5X[sds.72h.fem.5x > 2]

pm.72h.mal.5x <- percMethylation(DMS.meth.72h.mal.5X) #get percent methylation matrix
sds.72h.mal.5x <- matrixStats::rowSds(pm.72h.mal.5x) #calculate standard deviation of CpGs 
DMS.meth.72h.mal.5X <- DMS.meth.72h.mal.5X[sds.72h.mal.5x > 2]


#filter out SNPs
snp <- read.csv("../../BS-SNPer/shortterm_CT_SNP_edit.csv") #read in snps
snp.granges <- makeGRangesFromDataFrame(snp, ignore.strand = TRUE) #convert to granges 

DMS.meth.72h.5X <- DMS.meth.72h.5X[!as(DMS.meth.72h.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap
DMS.meth.72h.fem.5X <- DMS.meth.72h.fem.5X[!as(DMS.meth.72h.fem.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap
DMS.meth.72h.mal.5X <- DMS.meth.72h.mal.5X[!as(DMS.meth.72h.mal.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap

# Check number of CpGs 
DMS.meth.72h.5X
DMS.meth.72h.fem.5X
DMS.meth.72h.mal.5X

#enter covariates 
covariates.72h <- data.frame(tank=c("AC11","AC11","AC5","AC5","AC6","AC6",
                                    "C11","C11","C5","C5","C6","C6"), 
                             sex=c("F","M","F","M","F","M",
                                   "F","M","F","M","F","M"), 
                             stringsAsFactors = TRUE)

#calculate differential methylation
DMS.myDiff.72h.5X <- calculateDiffMeth(DMS.meth.72h.5X, covariates=covariates.72h, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMS_myDiff_72h")
DMS.myDiff.72h.fem.5X <- calculateDiffMeth(DMS.meth.72h.fem.5X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMS_myDiff_72h_fem")
DMS.myDiff.72h.mal.5X <- calculateDiffMeth(DMS.meth.72h.mal.5X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMS_myDiff_72h_mal")

#call significant methylation
DMS.diffMeth.72h.5X <- getMethylDiff(DMS.myDiff.72h.5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_72h")
DMS.diffMeth.72h.fem.5X <- getMethylDiff(DMS.myDiff.72h.fem.5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_72h_fem")
DMS.diffMeth.72h.mal.5X <- getMethylDiff(DMS.myDiff.72h.mal.5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_72h_mal")

#check number of significant DMS
DMS.diffMeth.72h.5X
DMS.diffMeth.72h.fem.5X
DMS.diffMeth.72h.mal.5X

# Get meth per chromosome
DMS.diffMethChr.72h.5X <- diffMethPerChr(DMS.myDiff.72h.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=15, save.db = TRUE, suffix = "DMS_chr_72h")
DMS.diffMethChr.72h.5X
DMS.diffMethChr.72h.fem.5X <- diffMethPerChr(DMS.myDiff.72h.fem.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=15, save.db = TRUE, suffix = "DMS_chr_72h_fem")
DMS.diffMethChr.72h.fem.5X
DMS.diffMethChr.72h.mal.5X <- diffMethPerChr(DMS.myDiff.72h.mal.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=15, save.db = TRUE, suffix = "DMS_chr_72h_mal")
DMS.diffMethChr.72h.mal.5X

## Save R objects ##
saveRDS(myobj.72h.subset, file = "./shortterm_myObj_72h_5X.RDS")
saveRDS(myobj.72h.fem.subset, file = "./shortterm_myObj_72h_fem_5X.RDS")
saveRDS(myobj.72h.mal.subset, file = "./shortterm_myObj_72h_mal_5X.RDS")

saveRDS(DMS.meth.72h.5X, file = "./shortterm_DMSmeth_72h_5X.RDS")
saveRDS(DMS.meth.72h.fem.5X, file = "./shortterm_DMSmeth_72h_fem_5X.RDS")
saveRDS(DMS.meth.72h.mal.5X, file = "./shortterm_DMSmeth_72h_mal_5X.RDS")

saveRDS(DMS.myDiff.72h.5X, file = "./shortterm_DMSmyDiff_72h_5X.RDS")
saveRDS(DMS.myDiff.72h.fem.5X, file = "./shortterm_DMSmyDiff_72h_fem_5X.RDS")
saveRDS(DMS.myDiff.72h.mal.5X, file = "./shortterm_DMSmyDiff_72h_mal_5X.RDS")

saveRDS(DMS.diffMeth.72h.5X, file = "./shortterm_DMSDiffMeth_72h_5X.RDS")
saveRDS(DMS.diffMeth.72h.fem.5X, file = "./shortterm_DMSDiffMeth_72h_fem_5X.RDS")
saveRDS(DMS.diffMeth.72h.mal.5X, file = "./shortterm_DMSDiffMeth_72h_mal_5X.RDS")

saveRDS(getData(DMS.diffMeth.72h.5X), file = "./shortterm_DMSDiffMeth_72h_5X_getData.RDS")
saveRDS(getData(DMS.diffMeth.72h.fem.5X), file = "./shortterm_DMSDiffMeth_72h_fem_5X_getData.RDS")
saveRDS(getData(DMS.diffMeth.72h.mal.5X), file = "./shortterm_DMSDiffMeth_72h_mal_5X_getData.RDS")

saveRDS(getData(DMS.myDiff.72h.5X), file = "./shortterm_DMSmyDiff_72h_5X_getData.RDS")
saveRDS(getData(DMS.myDiff.72h.fem.5X), file = "./shortterm_DMSmyDiff_72h_fem_5X_getData.RDS")
saveRDS(getData(DMS.myDiff.72h.mal.5X), file = "./shortterm_DMSmyDiff_72h_mal_5X_getData.RDS")
