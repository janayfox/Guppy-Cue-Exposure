#####################################################################################################################
### Goal: Read alignment files into methylKit and filter to create tabix files of filtered cytosine methylation, call DMS and DMRs
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

setwd("/scratch/janayfox/guppyWGBS/methylKit/dev/perc20_devSNP")

#Prepare function for getting DMS and DMRs 
run.DMS.DMR <- function(myobj,cov,chr.to.keep,unite.val,filename.myobj3x,filename.myobj5x,
                        file.name.meth.DMS,file.name.meth.DMS.data,
                        file.name.myDiff.DMS,file.name.myDiff.data,
                        file.name.diffMeth.DMS,filename.diffMeth.DMS.data,
                        filename.chr.DMS,
                        file.name.meth.DMR,file.name.meth.DMR.data,file.name.tile.meth,
                        file.name.myDiff.DMR,file.name.myDiff.DMR.data,
                        file.name.diffMeth.DMR,file.name.diffMeth.DMR.data,
                        filename.chr.DMR){

      #filter out sites in the 99.9th percentile of coverage (PCR bias) for 3x and 5X
      myobj.3X=filterByCoverage(myobj,lo.count=3,lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9)


      myobj.5X=filterByCoverage(myobj,lo.count=5,lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9)

      #normalize by median coverage
      norm.myobj.3X=normalizeCoverage(myobj.3X, method="median")
      norm.myobj.5X=normalizeCoverage(myobj.5X, method="median")

      #remove sex chr (LG12) and unplaced scaffolds
      myobj.subset.3X <- selectByOverlap(norm.myobj.3X, chr.to.keep)
      myobj.subset.5X <- selectByOverlap(norm.myobj.5X, chr.to.keep)

      ##Find DMS##
      #unite sites 
      DMS.meth.5X=unite(myobj.subset.5X, min.per.group = unite.val, destrand=FALSE, save.db = FALSE)

      #convert to non DB object 
      DMS.meth.5X <- as(DMS.meth.5X, "methylBase")

      #filter out low variation sites 
      pm.5X <- percMethylation(DMS.meth.5X) #get percent methylation matrix
      sds.5X <- matrixStats::rowSds(pm.5X, na.rm = TRUE) #calculate standard deviation of CpGs 
      DMS.meth.5X <- DMS.meth.5X[sds.5X > 2]

      #filter out SNPs
      snp <- read.csv("../../../BS-SNPer/dev_CT_SNP_edit.csv") #read in snps
      snp.granges <- makeGRangesFromDataFrame(snp, ignore.strand = TRUE) #convert to granges 
      DMS.meth.5X <- DMS.meth.5X[!as(DMS.meth.5X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap

      # Check number of CpGs 
      head(DMS.meth.5X)

      #calculate differential methylation
      DMS.myDiff.5X <- calculateDiffMeth(DMS.meth.5X, covariates=cov, mc.cores=2, test="Chisq", save.db = FALSE)

      #call significant methylation
      DMS.diffMeth.5X <- getMethylDiff(DMS.myDiff.5X, difference = 20, qvalue = 0.0125, save.db = FALSE)

      #check number of significant DMS
      head(DMS.diffMeth.5X)

      # Get meth per chromosome
      DMS.diffMeth.5X.chr <- diffMethPerChr(DMS.myDiff.5X, plot=FALSE, qvalue.cutoff=0.0125, meth.cutoff=20, save.db = FALSE)

      ##Find DMRs ##
      #unite sites 
      DMR.meth.3X=unite(myobj.subset.3X, min.per.group = unite.val, destrand=FALSE, save.db = FALSE)

      #convert to non DB object 
      DMR.meth.3X <- as(DMR.meth.3X, "methylBase")

      #filter out low variation sites 
      pm.3x <- percMethylation(DMR.meth.3X) #get percent methylation matrix
      sds.3x <- matrixStats::rowSds(pm.3x, na.rm = TRUE) #calculate standard deviation of CpGs 
      DMR.meth.3X <- DMR.meth.3X[sds.3x > 2]

      #filter out SNPs
      DMR.meth.3X <- DMR.meth.3X[!as(DMR.meth.3X, "GRanges") %over% snp.granges, ] #select CpGs that do not overlap

      # Check number of CpGs 
      head(DMR.meth.3X)

      #tile into 100 bp windows with min coverage 5X
      tile.meth.5X <- tileMethylCounts(DMR.meth.3X, win.size = 100, step.size = 100, cov.bases = 5)

      #check number of regions retained 
      head(tile.meth.5X)

      #calculate differential methylation 
      DMR.myDiff.5X <- calculateDiffMeth(tile.meth.5X, mc.cores=2, test="Chisq", covariates=cov, save.db = FALSE)
       
      #call significant methylation
      DMR.diffMeth.5X <- getMethylDiff(DMR.myDiff.5X, difference = 20, qvalue = 0.0125, save.db = FALSE)

      #check number of DMRs 
      head(DMR.diffMeth.5X)

      #get meth per chromosome 
      DMR.diffMeth.5X.chr <- diffMethPerChr(DMR.diffMeth.5X, plot = FALSE,qvalue.cutoff=0.0125, meth.cutoff=20, save.db = FALSE)

      # Save R objects ##
      saveRDS(myobj.subset.3X, file = filename.myobj3x)
      saveRDS(myobj.subset.5X, file = filename.myobj5x)

      saveRDS(DMS.meth.5X, file = file.name.meth.DMS)
      saveRDS(getData(DMS.meth.5X), file = file.name.meth.DMS.data)
      saveRDS(DMS.myDiff.5X, file = file.name.myDiff.DMS)
      saveRDS(getData(DMS.myDiff.5X), file = file.name.myDiff.data)
      saveRDS(DMS.diffMeth.5X, file = file.name.diffMeth.DMS)
      saveRDS(getData(DMS.diffMeth.5X), file = filename.diffMeth.DMS.data)
      saveRDS(DMS.diffMeth.5X.chr, file = filename.chr.DMS)

      
      saveRDS(DMR.meth.3X, file = file.name.meth.DMR)
      saveRDS(getData(DMR.meth.3X), file = file.name.meth.DMR.data)
      saveRDS(tile.meth.5X, file = file.name.tile.meth)
      saveRDS(DMR.myDiff.5X, file = file.name.myDiff.DMR)
      saveRDS(getData(DMR.myDiff.5X), file = file.name.myDiff.DMR.data)
      saveRDS(DMR.diffMeth.5X, file = file.name.diffMeth.DMR)
      saveRDS(getData(DMR.diffMeth.5X), file = file.name.diffMeth.DMR.data)
      saveRDS(DMR.diffMeth.5X.chr, file = filename.chr.DMR)

}

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

#read in tank covariate data
covariates.all <- data.frame(tank=c("AC2","AC2","AC2","AC3","AC3","AC3",
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
                              "F","F","M","M","F","F","F",
                              "M","M","M","F","F","F","F",
                              "F","M","M","F","F","F","F",
                              "F","M","M","M","F","F","F",
                              "F","M","M","M","F","F","F",
                              "M","M","M","M","M"),
                         stringsAsFactors = TRUE)

covariates.fem <- data.frame(tank=c("AC2","AC2","AC2","AC3","AC3",
                                "AC4","AC4","AC4","AC4","AC4",
                                "AC5","AC5","AC5","AC5",
                                "AC6","AC6","AC6",
                                "AC7","AC7","AC7",
                                "AC7","AC7","AC7",
                                "C2","C2","C3","C3",
                                "C3","C4","C4",
                                "C4","C4","C4","C5",
                                "C5","C5","C5","C5",
                                "C6","C6","C6","C6",
                                "C7","C7","C7"),
                         stringsAsFactors = TRUE)

covariates.mal <- data.frame(tank=c("AC3","AC3",
                                "AC5","AC5","AC5","AC5",
                                "AC6","AC6","AC6",
                                "AC7","AC7","AC7",
                                "C2","C2","C3","C3",
                                "C3","C4","C4",
                                "C5","C5","C5",
                                "C6","C6","C6",
                                "C7","C7","C7","C7",
                                "C7"),
                         stringsAsFactors = TRUE)

## Prepare tabix files
#create list of file locations
file.list.dev.all = list("../../../merged_sequences/dev/DAC2F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC2F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC2F6.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC3F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC3F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC3M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC3M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC4F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC4F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC4F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC4F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC4F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5M3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5M4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6M3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F6.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7M3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC2F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC2F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC2M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC2M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3M4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5M3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6M4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7M3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7M4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7M5.CpG_merged.cov")

file.list.dev.fem = list("../../../merged_sequences/dev/DAC2F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC2F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC2F6.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC3F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC3F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC4F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC4F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC4F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC4F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC4F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7F6.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC2F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC2F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5F5.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6F3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6F4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7F1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7F2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7F3.CpG_merged.cov")

file.list.dev.mal = list("../../../merged_sequences/dev/DAC3M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC3M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5M3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC5M4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC6M3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DAC7M3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC2M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC2M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC3M4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC4M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC5M3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC6M4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7M1.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7M2.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7M3.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7M4.CpG_merged.cov",
                     "../../../merged_sequences/dev/DC7M5.CpG_merged.cov")

#create tabix file
# myobj.all=methRead(file.list.dev.all,
#                sample.id=list("DAC2F4","DAC2F5","DAC2F6","DAC3F1","DAC3F2","DAC3M1",
#                               "DAC3M2","DAC4F1","DAC4F2","DAC4F3","DAC4F4","DAC4F5",
#                               "DAC5F1","DAC5F2","DAC5F4","DAC5F5","DAC5M1","DAC5M2",
#                               "DAC5M3","DAC5M4","DAC6F1","DAC6F2","DAC6F3","DAC6M1",
#                               "DAC6M2","DAC6M3", "DAC7F1","DAC7F2","DAC7F3","DAC7F4",
#                               "DAC7F5","DAC7F6","DAC7M1","DAC7M2","DAC7M3",
#                               "DC2F1","DC2F2", "DC2M1","DC2M2","DC3F1","DC3F2","DC3F3",
#                               "DC3M1","DC3M2","DC3M4","DC4F1","DC4F2","DC4F3","DC4F4",
#                               "DC4F5","DC4M1","DC4M2","DC5F1","DC5F2","DC5F3","DC5F4",
#                               "DC5F5","DC5M1","DC5M2","DC5M3","DC6F1","DC6F2","DC6F3",
#                               "DC6F4","DC6M1","DC6M2","DC6M4","DC7F1","DC7F2","DC7F3",
#                               "DC7M1","DC7M2","DC7M3","DC7M4","DC7M5"),
#                assembly="guppyWGBS_dev_final",
#                pipeline="bismarkCoverage",
#                treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
#                            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
#                context="CpG",
#               # dbtype = "tabix",
#                dbdir = "guppy_dev_DB",
#                mincov = 1
# )

# run.DMS.DMR(myobj.all,covariates.all,keep.chr.noXY,21L,
#                         "./DMS_res/myobj_all_3X.RDS","./DMS_res/myobj_all_5X.RDS",
#                         "./DMS_res/DMSmeth_all_5X.RDS","./DMS_res/meth_all_5X_data.RDS",
#                         "./DMS_res/DMSmydiff_all_5X.RDS","./DMS_res/DMSmydiff_all_5X_data.RDS",
#                         "./DMS_res/DMSdiffmeth_all_5X.RDS","./DMS_res/DMSdiffmeth_all_5X_data.RDS",
#                         "./DMS_res/DMSdiffmethchr_all_5X.RDS",
#                         "./DMR_res/DMRmeth_all_5X.RDS","./DMR_res/DMRmeth_all_5X_data.RDS",
#                         "./DMR_res/DMR_tile_meth_all_5X.RDS",
#                         "./DMR_res/DMRmydiff_all_5X.RDS","./DMR_res/DMRmydiff_all_5X_data.RDS",
#                         "./DMR_res/DMRdiffmeth_all_5X.RDS","./DMR_res/DMRdiffmeth_all_5X_data.RDS",
#                         "./DMR_res/DMRdiffmethchr_all_5X.RDS")

# rm(myobj.all)

# myobj.fem=methRead(file.list.dev.fem,
#                sample.id=list("DAC2F4","DAC2F5","DAC2F6","DAC3F1","DAC3F2",
#                               "DAC4F1","DAC4F2","DAC4F3","DAC4F4","DAC4F5",
#                               "DAC5F1","DAC5F2","DAC5F4","DAC5F5",
#                               "DAC6F1","DAC6F2","DAC6F3",
#                               "DAC7F1","DAC7F2","DAC7F3","DAC7F4",
#                               "DAC7F5","DAC7F6","DC2F1","DC2F2","DC3F1","DC3F2","DC3F3",
#                               "DC4F1","DC4F2","DC4F3","DC4F4",
#                               "DC4F5","DC5F1","DC5F2","DC5F3","DC5F4",
#                               "DC5F5","DC6F1","DC6F2","DC6F3",
#                               "DC6F4","DC7F1","DC7F2","DC7F3"),
#                assembly="guppyWGBS_dev_fem",
#                pipeline="bismarkCoverage",
#                treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
#                            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
#                context="CpG",
#               # dbtype = "tabix",
#                dbdir = "guppy_dev_DB_fem",
#                mincov = 1
# )

# run.DMS.DMR(myobj.fem,covariates.fem,keep.chr.allchr,13L,
#                         "./DMS_res/myobj_fem_3X.RDS","./DMS_res/myobj_fem_5X.RDS",
#                         "./DMS_res/DMSmeth_fem_5X.RDS","./DMS_res/meth_fem_5X_data.RDS",
#                         "./DMS_res/DMSmydiff_fem_5X.RDS","./DMS_res/DMSmydiff_fem_5X_data.RDS",
#                         "./DMS_res/DMSdiffmeth_fem_5X.RDS","./DMS_res/DMSdiffmeth_fem_5X_data.RDS",
#                         "./DMS_res/DMSdiffmethchr_fem_5X.RDS",
#                         "./DMR_res/DMRmeth_fem_5X.RDS","./DMR_res/DMRmeth_fem_5X_data.RDS",
#                         "./DMR_res/DMR_tile_meth_fem_5X.RDS",
#                         "./DMR_res/DMRmydiff_fem_5X.RDS","./DMR_res/DMRmydiff_fem_5X_data.RDS",
#                         "./DMR_res/DMRdiffmeth_fem_5X.RDS","./DMR_res/DMRdiffmeth_fem_5X_data.RDS",
#                         "./DMR_res/DMRdiffmethchr_fem_5X.RDS")

# rm(myobj.fem)

myobj.mal=methRead(file.list.dev.mal,
               sample.id=list("DAC3M1","DAC3M2","DAC5M1","DAC5M2",
                              "DAC5M3","DAC5M4","DAC6M1",
                              "DAC6M2","DAC6M3", "DAC7M1","DAC7M2","DAC7M3",
                              "DC2M1","DC2M2","DC3M1","DC3M2","DC3M4",
                              "DC4M1","DC4M2","DC5M1","DC5M2","DC5M3",
                              "DC6M1","DC6M2","DC6M4",
                              "DC7M1","DC7M2","DC7M3","DC7M4","DC7M5"),
               assembly="guppyWGBS_dev_mal",
               pipeline="bismarkCoverage",
               treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
               context="CpG",
              # dbtype = "tabix",
               dbdir = "guppy_dev_DB_mal",
               mincov = 1
)

run.DMS.DMR(myobj.mal,covariates.mal,keep.chr.allchr,7L,
                        "./DMS_res/myobj_mal_3X.RDS","./DMS_res/myobj_mal_5X.RDS",
                        "./DMS_res/DMSmeth_mal_5X.RDS","./DMS_res/meth_mal_5X_data.RDS",
                        "./DMS_res/DMSmydiff_mal_5X.RDS","./DMS_res/DMSmydiff_mal_5X_data.RDS",
                        "./DMS_res/DMSdiffmeth_mal_5X.RDS","./DMS_res/DMSdiffmeth_mal_5X_data.RDS",
                        "./DMS_res/DMSdiffmethchr_mal_5X.RDS",
                        "./DMR_res/DMRmeth_mal_5X.RDS","./DMR_res/DMRmeth_mal_5X_data.RDS",
                        "./DMR_res/DMR_tile_meth_mal_5X.RDS",
                        "./DMR_res/DMRmydiff_mal_5X.RDS","./DMR_res/DMRmydiff_mal_5X_data.RDS",
                        "./DMR_res/DMRdiffmeth_mal_5X.RDS","./DMR_res/DMRdiffmeth_mal_5X_data.RDS",
                        "./DMR_res/DMRdiffmethchr_mal_5X.RDS")
