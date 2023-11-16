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

setwd("/scratch/janayfox/guppyWGBS/methylKit/st/perc20/4h")

## Prepare tabix files
#create lists of file locations
file.list.4h = list("../../../../mergedCov/st/ST2AC13F.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2AC13M.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2AC1F.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2AC1M.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2AC7F.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2AC7M.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2C12F.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2C12M.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2C1F.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2C1M.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2C7F.CpG_merged.cov",
                    "../../../../mergedCov/st/ST2C7M.CpG_merged.cov")

file.list.4h.fem = list("../../../../mergedCov/st/ST2AC13F.CpG_merged.cov",
                        "../../../../mergedCov/st/ST2AC1F.CpG_merged.cov",
                        "../../../../mergedCov/st/ST2AC7F.CpG_merged.cov",
                        "../../../../mergedCov/st/ST2C12F.CpG_merged.cov",
                        "../../../../mergedCov/st/ST2C1F.CpG_merged.cov",
                        "../../../../mergedCov/st/ST2C7F.CpG_merged.cov")

file.list.4h.mal = list("../../../../mergedCov/st/ST2AC13M.CpG_merged.cov",
                        "../../../../mergedCov/st/ST2AC1M.CpG_merged.cov",
                        "../../../../mergedCov/st/ST2AC7M.CpG_merged.cov",
                        "../../../../mergedCov/st/ST2C12M.CpG_merged.cov",
                        "../../../../mergedCov/st/ST2C1M.CpG_merged.cov",
                        "../../../../mergedCov/st/ST2C7M.CpG_merged.cov")

#create tabix file
myobj.4h=methRead(file.list.4h,
                  sample.id=list("ST2AC13F","ST2AC13M","ST2AC1F","ST2AC1M","ST2AC7F","ST2AC7M",
                                 "ST2C12F","ST2C12M","ST2C1F","ST2C1M","ST2C7F","ST2C7M"),
                  assembly="guppyWGBS_shortterm",
                  pipeline="bismarkCoverage",
                  treatment=c(1,1,1,1,1,1,
                              0,0,0,0,0,0),
                  context="CpG",
                  mincov = 1,
                  dbtype = "tabix",
                  dbdir = "shortterm_4h_DB"
)

myobj.4h.fem=methRead(file.list.4h.fem,
                      sample.id=list("ST2AC13F","ST2AC1F","ST2AC7F",
                                     "ST2C12F","ST2C1F","ST2C7F"),
                      assembly="guppyWGBS_shortterm",
                      pipeline="bismarkCoverage",
                      treatment=c(1,1,1,
                                  0,0,0),
                      context="CpG",
                      mincov = 1,
                      dbtype = "tabix",
                      dbdir = "shortterm_4hF_DB"
)

myobj.4h.mal=methRead(file.list.4h.mal,
                      sample.id=list("ST2AC13M","ST2AC1M","ST2AC7M",
                                     "ST2C12M","ST2C1M","ST2C7M"),
                      assembly="guppyWGBS_shortterm",
                      pipeline="bismarkCoverage",
                      treatment=c(1,1,1,
                                  0,0,0),
                      context="CpG",
                      mincov = 1,
                      dbtype = "tabix",
                      dbdir = "shortterm_4hmal_DB"
)

#get coverage stats 
getCoverageStats(myobj.4h[[2]], both.strands = FALSE)
getCoverageStats(myobj.4h.fem[[2]], both.strands = FALSE)
getCoverageStats(myobj.4h.mal[[2]], both.strands = FALSE)

#filter out sites in the 99.9th percentile of coverage (PCR bias) 
myobj.4h.5X=filterByCoverage(myobj.4h,lo.count=5,lo.perc=NULL,
                               hi.count=NULL, hi.perc=99.9, suffix = "5X")
myobj.4h.fem.5X=filterByCoverage(myobj.4h.fem,lo.count=5,lo.perc=NULL,
                                   hi.count=NULL, hi.perc=99.9, suffix = "5X")
myobj.4h.mal.5X=filterByCoverage(myobj.4h.mal,lo.count=5,lo.perc=NULL,
                                   hi.count=NULL, hi.perc=99.9, suffix = "5X")

#normalize by median coverage
norm.myobj.4h.5X=normalizeCoverage(myobj.4h.5X, method="median")
norm.myobj.4h.fem.5X=normalizeCoverage(myobj.4h.fem.5X, method="median")
norm.myobj.4h.mal.5X=normalizeCoverage(myobj.4h.mal.5X, method="median")

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
myobj.4h.subset <- selectByOverlap(norm.myobj.4h.5X, keep.chr.noXY)
myobj.4h.fem.subset <- selectByOverlap(norm.myobj.4h.fem.5X, keep.chr.allchr)
myobj.4h.mal.subset <- selectByOverlap(norm.myobj.4h.mal.5X, keep.chr.allchr)

## Find DMS ##
#unite sites 
DMS.meth.4h.5X=unite(myobj.4h.subset, destrand=FALSE, min.per.group = 6L, save.db = TRUE, suffix = "DMS_unite_4h")
DMS.meth.4h.fem.5X=unite(myobj.4h.fem.subset, destrand=FALSE, min.per.group = 3L, save.db = TRUE, suffix = "DMS_unite_4h_fem")
DMS.meth.4h.mal.5X=unite(myobj.4h.mal.subset, destrand=FALSE, min.per.group = 3L, save.db = TRUE, suffix = "DMS_unite_4h_mal")

#convert to non DB object 
DMS.meth.4h.5X <- as(DMS.meth.4h.5X, "methylBase")
DMS.meth.4h.fem.5X <- as(DMS.meth.4h.fem.5X, "methylBase")
DMS.meth.4h.mal.5X <- as(DMS.meth.4h.mal.5X, "methylBase")

#filter out low variation sites 
pm.4h.5x <- percMethylation(DMS.meth.4h.5X) #get percent methylation matrix
sds.4h.5x <- matrixStats::rowSds(pm.4h.5x) #calculate standard deviation of CpGs 
DMS.meth.4h.5X <- DMS.meth.4h.5X[sds.4h.5x > 2]

pm.4h.fem.5x <- percMethylation(DMS.meth.4h.fem.5X) #get percent methylation matrix
sds.4h.fem.5x <- matrixStats::rowSds(pm.4h.fem.5x) #calculate standard deviation of CpGs 
DMS.meth.4h.fem.5X <- DMS.meth.4h.fem.5X[sds.4h.fem.5x > 2]

pm.4h.mal.5x <- percMethylation(DMS.meth.4h.mal.5X) #get percent methylation matrix
sds.4h.mal.5x <- matrixStats::rowSds(pm.4h.mal.5x) #calculate standard deviation of CpGs 
DMS.meth.4h.mal.5X <- DMS.meth.4h.mal.5X[sds.4h.mal.5x > 2]

#filter out SNPs
snp <- read.csv("../../../BS-SNPer/shortterm_CT_SNP_edit.csv") #read in snps
snp.granges <- makeGRangesFromDataFrame(snp, ignore.strand = TRUE) #convert to granges 

DMS.meth.4h.5X <- DMS.meth.4h.5X[!as(DMS.meth.4h.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap
DMS.meth.4h.fem.5X <- DMS.meth.4h.fem.5X[!as(DMS.meth.4h.fem.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap
DMS.meth.4h.mal.5X <- DMS.meth.4h.mal.5X[!as(DMS.meth.4h.mal.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap

#Check number of CpGs 
DMS.meth.4h.5X
DMS.meth.4h.fem.5X
DMS.meth.4h.mal.5X

#enter covariates 
covariates.4h <- data.frame(tank=c("AC13","AC13","AC1","AC1","AC7","AC7",
                                   "C12","C12","C1","C1","C7","C7"), 
                            sex=c("F","M","F","M","F","M",
                                  "F","M","F","M","F","M"), 
                            stringsAsFactors = TRUE)

#calculate differential methylation
DMS.myDiff.4h.5X <- calculateDiffMeth(DMS.meth.4h.5X, covariates=covariates.4h, mc.cores=2, test = "Chisq", save.db = TRUE, suffix = "DMS_myDiff_4h")
DMS.myDiff.4h.fem.5X <- calculateDiffMeth(DMS.meth.4h.fem.5X, mc.cores=2, test = "Chisq", save.db = TRUE, suffix = "DMS_myDiff_4h_fem")
DMS.myDiff.4h.mal.5X <- calculateDiffMeth(DMS.meth.4h.mal.5X, mc.cores=2, test = "Chisq", save.db = TRUE, suffix = "DMS_myDiff_4h_mal")

#call significant methylation
DMS.diffMeth.4h.5X <- getMethylDiff(DMS.myDiff.4h.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_4h")
DMS.diffMeth.4h.fem.5X <- getMethylDiff(DMS.myDiff.4h.fem.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_4h_fem")
DMS.diffMeth.4h.mal.5X <- getMethylDiff(DMS.myDiff.4h.mal.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_4h_mal")

#Check number of significant DMS
DMS.diffMeth.4h.5X
DMS.diffMeth.4h.fem.5X
DMS.diffMeth.4h.mal.5X

# Get meth per chromosome
DMS.diffMethChr.4h.5X <- diffMethPerChr(DMS.myDiff.4h.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = TRUE, suffix = "chr")
DMS.diffMethChr.4h.5X
DMS.diffMethChr.4h.fem.5X <- diffMethPerChr(DMS.myDiff.4h.fem.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = TRUE, suffix = "chr")
DMS.diffMethChr.4h.fem.5X
DMS.diffMethChr.4h.mal.5X <- diffMethPerChr(DMS.myDiff.4h.mal.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = TRUE, suffix = "chr")
DMS.diffMethChr.4h.mal.5X

## Save R objects ##
saveRDS(myobj.4h.subset, file = "./DMS_res/myobj_4h_5X.RDS")
saveRDS(myobj.4h.fem.subset, file = "./DMS_res/myobj_4h_fem_5X.RDS")
saveRDS(myobj.4h.mal.subset, file = "./DMS_res/myobj_4h_mal_5X.RDS")

saveRDS(DMS.meth.4h.5X, file = "./DMS_res/DMSmeth_4h_5X.RDS")
saveRDS(DMS.meth.4h.fem.5X, file = "./DMS_res/DMSmeth_4h_fem_5X.RDS")
saveRDS(DMS.meth.4h.mal.5X, file = "./DMS_res/DMSmeth_4h_mal_5X.RDS")

saveRDS(DMS.myDiff.4h.5X, file = "./DMS_res/DMSmydiff_4h_5X.RDS")
saveRDS(DMS.myDiff.4h.fem.5X, file = "./DMS_res/DMSmydiff_4h_fem_5X.RDS")
saveRDS(DMS.myDiff.4h.mal.5X, file = "./DMS_res/DMSmydiff_4h_mal_5X.RDS")

saveRDS(DMS.diffMeth.4h.5X, file = "./DMS_res/DMSdiffmeth_4h_5X.RDS")
saveRDS(DMS.diffMeth.4h.fem.5X, file = "./DMS_res/DMSdiffmeth_4h_fem_5X.RDS")
saveRDS(DMS.diffMeth.4h.mal.5X, file = "./DMS_res/DMSdiffmeth_4h_mal_5X.RDS")

saveRDS(getData(DMS.meth.4h.5X), file = "./DMS_res/meth_4h_5X_data.RDS")
saveRDS(getData(DMS.meth.4h.fem.5X), file = "./DMS_res/meth_4h_fem_5X_data.RDS")
saveRDS(getData(DMS.meth.4h.mal.5X), file = "./DMS_res/meth_4h_mal_5X_data.RDS")

saveRDS(getData(DMS.diffMeth.4h.5X), file = "./DMS_res/DMSdiffmeth_4h_5X_data.RDS")
saveRDS(getData(DMS.diffMeth.4h.fem.5X), file = "./DMS_res/DMSdiffmeth_4h_fem_5X_data.RDS")
saveRDS(getData(DMS.diffMeth.4h.mal.5X), file = "./DMS_res/DMSdiffmeth_4h_mal_5X_data.RDS")

saveRDS(getData(DMS.myDiff.4h.5X), file = "./DMS_res/DMSmydiff_4h_5X_data.RDS")
saveRDS(getData(DMS.myDiff.4h.fem.5X), file = "./DMS_res/DMSmydiff_4h_fem_5X_data.RDS")
saveRDS(getData(DMS.myDiff.4h.mal.5X), file = "./DMS_res/DMSmydiff_4h_mal_5X_data.RDS")
