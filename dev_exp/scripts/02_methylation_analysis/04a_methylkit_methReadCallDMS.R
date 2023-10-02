#####################################################################################################################
### Goal: Read alignment files into methylKit and filter to create tabix files of filtered cytosine methylation, call DMS
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
library("data.table", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

setwd("/scratch/janayfox/guppyWGBS/methylKit/dev")

## Prepare tabix files
#create list of file locations
file.list.dev = list("../../mergedCov/dev/DAC2F4.CpG_merged.cov",
                     "../../mergedCov/dev/DAC2F5.CpG_merged.cov",
                     "../../mergedCov/dev/DAC2F6.CpG_merged.cov",
                     "../../mergedCov/dev/DAC3F1.CpG_merged.cov",
                     "../../mergedCov/dev/DAC3F2.CpG_merged.cov",
                     "../../mergedCov/dev/DAC3M1.CpG_merged.cov",
                     "../../mergedCov/dev/DAC3M2.CpG_merged.cov",
                     "../../mergedCov/dev/DAC4F1.CpG_merged.cov",
                     "../../mergedCov/dev/DAC4F2.CpG_merged.cov",
                     "../../mergedCov/dev/DAC4F3.CpG_merged.cov",
                     "../../mergedCov/dev/DAC4F4.CpG_merged.cov",
                     "../../mergedCov/dev/DAC4F5.CpG_merged.cov",
                     "../../mergedCov/dev/DAC5F1.CpG_merged.cov",
                     "../../mergedCov/dev/DAC5F2.CpG_merged.cov",
                     "../../mergedCov/dev/DAC5F4.CpG_merged.cov",
                     "../../mergedCov/dev/DAC5F5.CpG_merged.cov",
                     "../../mergedCov/dev/DAC5M1.CpG_merged.cov",
                     "../../mergedCov/dev/DAC5M2.CpG_merged.cov",
                     "../../mergedCov/dev/DAC5M3.CpG_merged.cov",
                     "../../mergedCov/dev/DAC5M4.CpG_merged.cov",
                     "../../mergedCov/dev/DAC6F1.CpG_merged.cov",
                     "../../mergedCov/dev/DAC6F2.CpG_merged.cov",
                     "../../mergedCov/dev/DAC6F3.CpG_merged.cov",
                     "../../mergedCov/dev/DAC6M1.CpG_merged.cov",
                     "../../mergedCov/dev/DAC6M2.CpG_merged.cov",
                     "../../mergedCov/dev/DAC6M3.CpG_merged.cov",
                     "../../mergedCov/dev/DAC7F1.CpG_merged.cov",
                     "../../mergedCov/dev/DAC7F2.CpG_merged.cov",
                     "../../mergedCov/dev/DAC7F3.CpG_merged.cov",
                     "../../mergedCov/dev/DAC7F4.CpG_merged.cov",
                     "../../mergedCov/dev/DAC7F5.CpG_merged.cov",
                     "../../mergedCov/dev/DAC7F6.CpG_merged.cov",
                     "../../mergedCov/dev/DAC7M1.CpG_merged.cov",
                     "../../mergedCov/dev/DAC7M2.CpG_merged.cov",
                     "../../mergedCov/dev/DAC7M3.CpG_merged.cov",
                     "../../mergedCov/dev/DC2F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC2F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC2M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC2M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC3F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC3M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC3M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC3M4.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F4.CpG_merged.cov",
                     "../../mergedCov/dev/DC4F5.CpG_merged.cov",
                     "../../mergedCov/dev/DC4M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC4M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F4.CpG_merged.cov",
                     "../../mergedCov/dev/DC5F5.CpG_merged.cov",
                     "../../mergedCov/dev/DC5M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC5M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC5M3.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC6F4.CpG_merged.cov",
                     "../../mergedCov/dev/DC6M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC6M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC6M4.CpG_merged.cov",
                     "../../mergedCov/dev/DC7F1.CpG_merged.cov",
                     "../../mergedCov/dev/DC7F2.CpG_merged.cov",
                     "../../mergedCov/dev/DC7F3.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M1.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M2.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M3.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M4.CpG_merged.cov",
                     "../../mergedCov/dev/DC7M5.CpG_merged.cov")

#create tabix file
myobj=methRead(file.list.dev,
               sample.id=list("DAC2F4","DAC2F5","DAC2F6","DAC3F1","DAC3F2","DAC3M1",
                              "DAC3M2","DAC4F1","DAC4F2","DAC4F3","DAC4F4","DAC4F5",
                              "DAC5F1","DAC5F2","DAC5F4","DAC5F5","DAC5M1","DAC5M2",
                              "DAC5M3","DAC5M4","DAC6F1","DAC6F2","DAC6F3","DAC6M1",
                              "DAC6M2","DAC6M3", "DAC7F1","DAC7F2","DAC7F3","DAC7F4",
                              "DAC7F5","DAC7F6","DAC7M1","DAC7M2","DAC7M3",
                              "DC2F1","DC2F2", "DC2M1","DC2M2","DC3F1","DC3F2","DC3F3",
                              "DC3M1","DC3M2","DC3M4","DC4F1","DC4F2","DC4F3","DC4F4",
                              "DC4F5","DC4M1","DC4M2","DC5F1","DC5F2","DC5F3","DC5F4",
                              "DC5F5","DC5M1","DC5M2","DC5M3","DC6F1","DC6F2","DC6F3",
                              "DC6F4","DC6M1","DC6M2","DC6M4","DC7F1","DC7F2","DC7F3",
                              "DC7M1","DC7M2","DC7M3","DC7M4","DC7M5"),
               assembly="guppyWGBS_dev_final",
               pipeline="bismarkCoverage",
               treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
               context="CpG",
               dbtype = "tabix",
               dbdir = "guppy_dev_DB_final",
               mincov = 1
)

#get coverage stats 
getCoverageStats(myobj[[2]], both.strands = FALSE)

#filter out sites in the 99.9th percentile of coverage (PCR bias) and 5x coverage
myobj.5X=filterByCoverage(myobj,lo.count=5,lo.perc=NULL,
                          hi.count=NULL, hi.perc=99.9, suffix = "5X_merged")

#normalize by median coverage
norm.myobj.5X=normalizeCoverage(myobj.5X, method="median")

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

myobj5X.subset <- selectByOverlap(norm.myobj.5X, keep.chr.noXY)

## Find DMSs ## 
# Unite methylation calls
DMSmeth5X_21L <- unite(myobj5X.subset, min.per.group=21L, destrand=FALSE, save.db = TRUE, suffix = "DMSunite5X_21L")

#filter out low variation sites 
pm <- percMethylation(DMSmeth5X_21L) #get percent methylation matrix
sds <- matrixStats::rowSds(pm) #calculate standard deviation of CpGs 
DMSmeth5X_21L <- DMSmeth5X_21L[sds > 2]

#filter out SNPs
snp <- read.csv("../../BS-SNPer/dev_CT_SNP_edit.csv") #read in snps
snp.granges <- makeGRangesFromDataFrame(snp, ignore.strand = TRUE) #convert to granges 
DMSmeth5X_21L <- DMSmeth5X_21L[!as(DMSmeth5X_21L, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap

#check number of cpgs
nrow(DMSmeth5X_21L)

#read in tank covariate data
covariates <- data.frame(tank=c("AC2","AC2","AC2","AC3","AC3","AC3",
                                "AC3","AC4","AC4","AC4","AC4","AC4",
                                "AC5","AC5","AC5","AC5","AC5",
                                "AC5","AC5","AC5","AC6","AC6","AC6",
                                "AC6","AC6","AC6","AC7","AC7","AC7",
                                "AC7","AC7","AC7","AC7","AC7","AC7",
                                "C2","C2","C2","C2","C3","C3",
                                "C3","C3","C3","C3","C4","C4",
                                "C4","C4","C4","C4","C4","C5",
                                "C5","C5","C5","C5","C5","C5",
                                "C5","C6","C6","C6","C6","C6",
                                "C6","C6","C7","C7","C7","C7",
                                "C7","C7","C7","C7"),
                        sex=c("F","F","F","F","F","M",
                              "M","F","F","F","F","F",
                              "F","F","F","F","M","M",
                              "M","M","F","F","F","M",
                              "M","M", "F","F","F","F",
                              "F","F","M","M","M",
                              "F","F", "M","M","F","F","F",
                              "M","M","M","F","F","F","F",
                              "F","M","M","F","F","F","F",
                              "F","M","M","M","F","F","F",
                              "F","M","M","M","F","F","F",
                              "M","M","M","M","M"),
                         stringsAsFactors = TRUE)

# Calculate differential methylation
DMSmyDiff5X_21L <- calculateDiffMeth(DMSmeth5X_21L, mc.cores=2, covariates=covariates, test="Chisq", save.db = TRUE, suffix = "myDiff_merged")

# Call significant methylation
DMSdiffMeth5X_21L <- getMethylDiff(DMSmyDiff5X_21L, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "diffMeth_merged")

# Check number of significant DMS
DMSdiffMeth5X_21L

# Get meth per chromosome
DMSdiffMethChr5X_21L <- diffMethPerChr(DMSmyDiff5X_21L, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=15, save.db = TRUE, suffix = "chr_merged")
DMSdiffMethChr5X_21L

## Save R objects ##
saveRDS(myobj5X.subset, file = "./myObj5X.RDS")
saveRDS(DMSmeth5X_21L, file = "./DMSmeth5x_25L.RDS")
saveRDS(DMSmyDiff5X_21L, file = "./DMSmyDiff5X_25L.RDS")
saveRDS(DMSdiffMeth5X_21L, file = "./DMSdiffmeth5X_25L.RDS")
saveRDS(DMSdiffMethChr5X_21L, file = "./DMSdiffMethChr5x_25L.RDS")
saveRDS(getData(DMSdiffMeth5X_21L), file = "./DMSdiffMeth5X_25Ldata.RDS")
saveRDS(getData(DMSmyDiff5X_25L), file = "./DMSmyDiff5X_21Ldata.RDS")


