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

setwd("/scratch/janayfox/guppyWGBS/methylKit/st/perc20/all/")

## Prepare tabix files
#create lists of file locations
file.list.all = list("../../../../mergedCov/st/ST2AC10F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC10M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC11F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC11M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC13F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC13M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC14F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC14M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC15F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC15M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC16F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC16M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC1F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC1M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC2F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC2M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC3F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC3M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC4F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC4M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC5F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC5M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC6F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC6M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC7F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC7M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC8F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC8M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC9F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2AC9M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C10F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C10M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C11F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C11M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C12F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C12M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C13F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C13M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C14F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C14M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C15F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C15M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C1F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C1M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C2F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C2M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C3F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C3M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C4F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C4M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C5F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C5M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C6F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C6M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C7F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C7M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C8F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C8M.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C9F.CpG_merged.cov",
                     "../../../../mergedCov/st/ST2C9M.CpG_merged.cov")

file.list.all.fem = list("../../../../mergedCov/st/ST2AC10F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC11F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC13F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC14F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC15F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC16F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC1F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC2F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC3F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC4F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC5F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC6F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC7F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC8F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC9F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C10F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C11F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C12F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C13F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C14F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C15F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C1F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C2F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C3F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C4F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C5F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C6F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C7F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C8F.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C9F.CpG_merged.cov")

file.list.all.mal = list("../../../../mergedCov/st/ST2AC10M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC11M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC13M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC14M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC15M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC16M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC1M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC2M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC3M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC4M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC5M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC6M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC7M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC8M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2AC9M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C10M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C11M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C12M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C13M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C14M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C15M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C1M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C2M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C3M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C4M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C5M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C6M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C7M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C8M.CpG_merged.cov",
                         "../../../../mergedCov/st/ST2C9M.CpG_merged.cov")

#create tabix file
myobj.all=methRead(file.list.all,
                   sample.id=list("ST2AC10F","ST2AC10M","ST2AC11F","ST2AC11M","ST2AC13F",
                                  "ST2AC13M","ST2AC14F","ST2AC14M","ST2AC15F","ST2AC15M",
                                  "ST2AC16F","ST2AC16M","ST2AC1F","ST2AC1M","ST2AC2F",
                                  "ST2AC2M","ST2AC3F","ST2AC3M","ST2AC4F","ST2AC4M",
                                  "ST2AC5F","ST2AC5M","ST2AC6F","ST2AC6M","ST2AC7F",
                                  "ST2AC7M","ST2AC8F","ST2AC8M","ST2AC9F","ST2AC9M",
                                  "ST2C10F","ST2C10M","ST2C11F","ST2C11M","ST2C12F",
                                  "ST2C12M","ST2C13F","ST2C13M","ST2C14F","ST2C14M",
                                  "ST2C15F","ST2C15M","ST2C1F","ST2C1M","ST2C2F",
                                  "ST2C2M","ST2C3F","ST2C3M","ST2C4F","ST2C4M",
                                  "ST2C5F","ST2C5M","ST2C6F","ST2C6M","ST2C7F",
                                  "ST2C7M","ST2C8F","ST2C8M","ST2C9F", "ST2C9M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCoverage",
                   treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                   context="CpG",
                   mincov = 1,
                   dbtype = "tabix",
                   dbdir = "shortterm_all_DB"
)

myobj.all.fem=methRead(file.list.all.fem,
                       sample.id=list("ST2AC10F","ST2AC11F","ST2AC13F","ST2AC14F","ST2AC15F",
                                      "ST2AC16F","ST2AC1F","ST2AC2F","ST2AC3F","ST2AC4F",
                                      "ST2AC5F","ST2AC6F","ST2AC7F","ST2AC8F","ST2AC9F",
                                      "ST2C10F","ST2C11F","ST2C12F","ST2C13F","ST2C14F",
                                      "ST2C15F","ST2C1F","ST2C2F","ST2C3F","ST2C4F",
                                      "ST2C5F","ST2C6F","ST2C7F","ST2C8F","ST2C9F"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_allfem_DB"
)

myobj.all.mal=methRead(file.list.all.mal,
                       sample.id=list("ST2AC10M","ST2AC11M","ST2AC13M","ST2AC14M","ST2AC15M",
                                      "ST2AC16M","ST2AC1M","ST2AC2M","ST2AC3M","ST2AC4M",
                                      "ST2AC5M","ST2AC6M","ST2AC7M","ST2AC8M","ST2AC9M",
                                      "ST2C10M","ST2C11M","ST2C12M","ST2C13M","ST2C14M",
                                      "ST2C15M","ST2C1M","ST2C2M","ST2C3M","ST2C4M",
                                      "ST2C5M","ST2C6M","ST2C7M","ST2C8M", "ST2C9M"),
                       assembly="guppyWGBS_shortterm",
                       pipeline="bismarkCoverage",
                       treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                       context="CpG",
                       mincov = 1,
                       dbtype = "tabix",
                       dbdir = "shortterm_allmal_DB"
)

#get coverage stats 
getCoverageStats(myobj.all[[2]], both.strands = FALSE)
getCoverageStats(myobj.all.fem[[2]], both.strands = FALSE)
getCoverageStats(myobj.all.mal[[2]], both.strands = FALSE)

#filter out sites in the 99.9th percentile of coverage (PCR bias) 
myobj.all.5X=filterByCoverage(myobj.all,lo.count=5,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj.all.fem.5X=filterByCoverage(myobj.all.fem,lo.count=5,lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9, suffix = "filt")
myobj.all.mal.5X=filterByCoverage(myobj.all.mal,lo.count=5,lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9, suffix = "filt")

#normalize by median coverage
norm.myobj.all.5X=normalizeCoverage(myobj.all.5X, method="median")
norm.myobj.all.fem.5X=normalizeCoverage(myobj.all.fem.5X, method="median")
norm.myobj.all.mal.5X=normalizeCoverage(myobj.all.mal.5X, method="median")

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

#remove sex chr (LG12) and unplaced scaffolds
myobj.all.subset <- selectByOverlap(norm.myobj.all.5X, keep.chr.noXY)
myobj.all.fem.subset <- selectByOverlap(norm.myobj.all.fem.5X, keep.chr.allchr)
myobj.all.mal.subset <- selectByOverlap(norm.myobj.all.mal.5X, keep.chr.allchr)

##Find DMS##
#unite sites 
DMS.meth.all.5X=unite(myobj.all.subset, min.per.group = 30L, destrand=FALSE, save.db = TRUE, suffix = "DMS_unite_all")
DMS.meth.all.fem.5X=unite(myobj.all.fem.subset, min.per.group = 15L, destrand=FALSE, save.db = TRUE, suffix = "DMS_unite_all_fem")
DMS.meth.all.mal.5X=unite(myobj.all.mal.subset, min.per.group = 15L, destrand=FALSE, save.db = TRUE, suffix = "DMS_unite_all_mal")

#convert to non DB object 
DMS.meth.all.5X <- as(DMS.meth.all.5X, "methylBase")
DMS.meth.all.fem.5X <- as(DMS.meth.all.fem.5X, "methylBase")
DMS.meth.all.mal.5X <- as(DMS.meth.all.mal.5X, "methylBase")

#filter out low variation sites 
pm.all.5x <- percMethylation(DMS.meth.all.5X) #get percent methylation matrix
sds.all.5x <- matrixStats::rowSds(pm.all.5x) #calculate standard deviation of CpGs 
DMS.meth.all.5X <- DMS.meth.all.5X[sds.all.5x > 2]

pm.all.fem.5x <- percMethylation(DMS.meth.all.fem.5X) #get percent methylation matrix
sds.all.fem.5x <- matrixStats::rowSds(pm.all.fem.5x) #calculate standard deviation of CpGs 
DMS.meth.all.fem.5X <- DMS.meth.all.fem.5X[sds.all.fem.5x > 2]

pm.all.mal.5x <- percMethylation(DMS.meth.all.mal.5X) #get percent methylation matrix
sds.all.mal.5x <- matrixStats::rowSds(pm.all.mal.5x) #calculate standard deviation of CpGs 
DMS.meth.all.mal.5X <- DMS.meth.all.mal.5X[sds.all.mal.5x > 2]

#filter out SNPs
snp <- read.csv("../../../../BS-SNPer/shortterm_CT_SNP_edit.csv") #read in snps
snp.granges <- makeGRangesFromDataFrame(snp, ignore.strand = TRUE) #convert to granges 

DMS.meth.all.5X <- DMS.meth.all.5X[!as(DMS.meth.all.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap
DMS.meth.all.fem.5X <- DMS.meth.all.fem.5X[!as(DMS.meth.all.fem.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap
DMS.meth.all.mal.5X <- DMS.meth.all.mal.5X[!as(DMS.meth.all.mal.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap

# Check number of CpGs 
DMS.meth.all.5X
DMS.meth.all.fem.5X
DMS.meth.all.mal.5X

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
                            sex=c("F","M","F","M","F",
                                  "M","F","M","F","M",
                                  "F","M","F","M","F",
                                  "M","F","M","F","M",
                                  "F","M","F","M","F",
                                  "M","F","M","F","M",
                                  "F","M","F","M","F",
                                  "M","F","M","F","M",
                                  "F","M","F","M","F",
                                  "M","F","M","F","M",
                                  "F","M","F","M","F",
                                  "M","F","M","F", "M"), 
                            stringsAsFactors = TRUE)

#calculate differential methylation
DMS.myDiff.all.5X <- calculateDiffMeth(DMS.meth.all.5X, covariates=covariates.all, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMS_myDiff_all")
DMS.myDiff.all.fem.5X <- calculateDiffMeth(DMS.meth.all.fem.5X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMS_myDiff_all_fem")
DMS.myDiff.all.mal.5X <- calculateDiffMeth(DMS.meth.all.mal.5X, mc.cores=2, test="Chisq", save.db = TRUE, suffix = "DMS_myDiff_all_mal")

#call significant methylation
DMS.diffMeth.all.5X <- getMethylDiff(DMS.myDiff.all.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_all")
DMS.diffMeth.all.fem.5X <- getMethylDiff(DMS.myDiff.all.fem.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_all_fem")
DMS.diffMeth.all.mal.5X <- getMethylDiff(DMS.myDiff.all.mal.5X, difference = 20, qvalue = 0.0125, save.db = TRUE, suffix = "DMS_diffMeth_all_mal")

#check number of significant DMS
DMS.diffMeth.all.5X
DMS.diffMeth.all.fem.5X
DMS.diffMeth.all.mal.5X

# Get meth per chromosome
DMS.diffMethChr.all.5X <- diffMethPerChr(DMS.myDiff.all.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = TRUE, suffix = "DMS_chr_all")
DMS.diffMethChr.all.5X
DMS.diffMethChr.all.fem.5X <- diffMethPerChr(DMS.myDiff.all.fem.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = TRUE, suffix = "DMS_chr_all_fem")
DMS.diffMethChr.all.fem.5X
DMS.diffMethChr.all.mal.5X <- diffMethPerChr(DMS.myDiff.all.mal.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = TRUE, suffix = "DMS_chr_all_mal")
DMS.diffMethChr.all.mal.5X

## Save R objects ##
saveRDS(myobj.all.subset, file = "./DMS_res/myobj_all_5X.RDS")
saveRDS(myobj.all.fem.subset, file = "./DMS_res/myobj_all_fem_5X.RDS")
saveRDS(myobj.all.mal.subset, file = "./DMS_res/myobj_all_mal_5X.RDS")

saveRDS(DMS.meth.all.5X, file = "./DMS_res/DMSmeth_all_5X.RDS")
saveRDS(DMS.meth.all.fem.5X, file = "./DMS_res/DMSmeth_all_fem_5X.RDS")
saveRDS(DMS.meth.all.mal.5X, file = "./DMS_res/DMSmeth_all_mal_5X.RDS")

saveRDS(DMS.myDiff.all.5X, file = "./DMS_res/DMSmydiff_all_5X.RDS")
saveRDS(DMS.myDiff.all.fem.5X, file = "./DMS_res/DMSmydiff_all_fem_5X.RDS")
saveRDS(DMS.myDiff.all.mal.5X, file = "./DMS_res/DMSmydiff_all_mal_5X.RDS")

saveRDS(DMS.diffMeth.all.5X, file = "./DMS_res/DMSdiffmeth_all_5X.RDS")
saveRDS(DMS.diffMeth.all.fem.5X, file = "./DMS_res/DMSdiffmeth_all_fem_5X.RDS")
saveRDS(DMS.diffMeth.all.mal.5X, file = "./DMS_res/DMSdiffmeth_all_mal_5X.RDS")

saveRDS(getData(DMS.meth.all.5X), file = "./DMS_res/meth_all_5X_data.RDS")
saveRDS(getData(DMS.meth.all.fem.5X), file = "./DMS_res/meth_all_fem_5X_data.RDS")
saveRDS(getData(DMS.meth.all.mal.5X), file = "./DMS_res/meth_all_mal_5X_data.RDS")

saveRDS(getData(DMS.diffMeth.all.5X), file = "./DMS_res/DMSdiffmeth_all_5X_data.RDS")
saveRDS(getData(DMS.diffMeth.all.fem.5X), file = "./DMS_res/DMSdiffmeth_all_fem_5X_data.RDS")
saveRDS(getData(DMS.diffMeth.all.mal.5X), file = "./DMS_res/DMSdiffmeth_all_mal_5X_data.RDS")

saveRDS(getData(DMS.myDiff.all.5X), file = "./DMS_res/DMSmydiff_all_5X_data.RDS")
saveRDS(getData(DMS.myDiff.all.fem.5X), file = "./DMS_res/DMSmydiff_all_fem_5X_data.RDS")
saveRDS(getData(DMS.myDiff.all.mal.5X), file = "./DMS_res/DMSmydiff_all_mal_5X_data.RDS")

