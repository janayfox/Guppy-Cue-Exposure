###########################################################################################
### Goal: Remove sex chr and unplaced scaffolds and methylation calling for DMSs and DMRs
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

setwd("/scratch/janayfox/guppyWGBS/")

#read in R obj 
myobj10X <- readRDS("./myObj10X.RDS")
myobj5X <- readRDS("./myObj5X.RDS")

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
                         stringsAsFactors = TRUE)

## Remove sex chromosomes ##
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

myobj5X.subset <- selectByOverlap(myobj5X, keep.chr.noXY)

#check how many CpGs I have now 
myobj5X.subset

## Find DMSs ## 
# Unite methylation calls for 70% of fish
#DMSmeth10X <- unite(myobj.10X, min.per.group=20L, destrand=TRUE, save.db = TRUE, suffix = "DMSunite10X")
DMSmeth5X <- unite(myobj5X.subset, min.per.group=25L, destrand=TRUE, save.db = TRUE, suffix = "DMSunite5X")

DMSmeth5X <- unite(myobj5X, min.per.group=25L, destrand=TRUE, save.db = TRUE, suffix = "DMSunite5X")


# Check number of CpGs 
#DMSmeth10X
DMSmeth5X

# Calculate differential methylation
#DMSmyDiff10X <- calculateDiffMeth(DMSmeth10X, mc.cores=2, covariates=covariates, overdispersion="MN", test="Chisq", save.db = TRUE, suffix = "myDiff")
DMSmyDiff5X <- calculateDiffMeth(DMSmeth5X, mc.cores=2, covariates=covariates, overdispersion="MN", test="Chisq", save.db = TRUE, suffix = "myDiff")

# Call significant methylation
#DMSdiffMeth10X <- getMethylDiff(DMSmyDiff10X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "diffMeth")
DMSdiffMeth5X <- getMethylDiff(DMSmyDiff5X, difference = 15, qvalue = 0.0125, save.db = TRUE, suffix = "diffMeth")

# Check number of significant DMS
#DMSdiffMeth10X
DMSdiffMeth5X

# Get meth per chromosome
#DMSdiffMethChr10X <- diffMethPerChr(DMSmyDiff10X,plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "chr")
#DMSdiffMethChr10X
DMSdiffMethChr5X <- diffMethPerChr(DMSmyDiff5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=15, save.db =  TRUE, suffix = "chr")
DMSdiffMethChr5X

## Save R objects ##
#saveRDS(DMSmeth10X, file = "./DMSmeth10X.RDS")
saveRDS(DMSmeth5X, file = "./DMSmeth5X.RDS")

#saveRDS(DMSmyDiff10X, file = "./DMSmyDiff10X.RDS")
saveRDS(DMSmyDiff5X, file = "./DMSmyDiff5X.RDS")

#saveRDS(DMSdiffMeth10X, file = "./DMSdiffMeth10X.RDS")
saveRDS(DMSdiffMeth5X, file = "./DMSdiffMeth5X.RDS")

#saveRDS(DMSdiffMethChr10X, file = "./DMSdiffMethChr10X.RDS")
saveRDS(DMSdiffMethChr5X, file = "./DMSdiffMethChr5X.RDS")

## Get data ## 
#getData(DMSmeth10X) %>% saveRDS(file = "./DMSmeth10Xdata.RDS")
#getData(DMSmyDiff10X) %>% saveRDS(file = "./DMSmyDiff10Xdata.RDS")
#getData(DMSdiffMeth10X) %>% saveRDS(file = "./DMSdiffMeth10Xdata.RDS")

saveRDS(getData(DMSdiffMeth5X), file = "./DMSdiffMeth5Xdata.RDS")
saveRDS(getData(DMSmyDiff5X), file = "./DMSmyDiff5Xdata.RDS")
saveRDS(getData(DMSmeth5X), file = "./DMSmeth5Xdata.RDS")
