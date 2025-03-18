##########################################################
### Goal: Do permutation test on DMS/DMR analysis 
### Author: Janay Fox
### R script
###########################################################

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

setwd("/scratch/janayfox/guppyWGBS_final/methylKit/st/perc20_permutation")

#Prepare function for getting DMR 
run.DMR <- function(myobj,cov,chr.to.keep,unite.val){
  
  #filter out sites in the 99.9th percentile of coverage (PCR bias) for 3x and 5X
  myobj.3X=filterByCoverage(myobj,lo.count=3,lo.perc=NULL,
                            hi.count=NULL, hi.perc=99.9)
  
  #normalize by median coverage
  norm.myobj.3X=normalizeCoverage(myobj.3X, method="median")
  
  #remove sex chr (LG12) and unplaced scaffolds
  myobj.subset.3X <- selectByOverlap(norm.myobj.3X, chr.to.keep)
  
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
  snp <- read.csv("../../../BS-SNPer/shortterm_CT_SNP_edit.csv") #read in snps
  snp.granges <- makeGRangesFromDataFrame(snp, ignore.strand = TRUE) #convert to granges   
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
  num.DMR <- nrow(DMR.diffMeth.5X)
  
  #return number of DMRs 
  return(num.DMR)
  
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

#create lists of file locations
file.list.all.fem = list("../../../merged_sequences/st/ST2AC10F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC11F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC13F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC14F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC15F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC16F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC1F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC2F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC3F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC4F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC5F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC6F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC7F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC8F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC9F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C10F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C11F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C12F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C13F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C14F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C15F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C1F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C2F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C3F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C4F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C5F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C6F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C7F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C8F.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C9F.CpG_merged.cov")

file.list.all.mal = list("../../../merged_sequences/st/ST2AC10M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC11M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC13M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC14M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC15M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC16M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC1M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC2M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC3M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC4M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC5M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC6M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC7M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC8M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2AC9M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C10M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C11M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C12M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C13M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C14M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C15M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C1M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C2M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C3M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C4M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C5M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C6M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C7M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C8M.CpG_merged.cov",
                         "../../../merged_sequences/st/ST2C9M.CpG_merged.cov")

#set up permutation test
set.seed(113)
nsim <- 100
num.DMR.fem <- numeric(nsim)
num.DMR.mal <- numeric(nsim)

#treatment vector
original.treatment.fem <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

original.treatment.mal <-  c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#run permutation test on females 
for (i in 1:nsim){
  
  #randomly shuffle the treatment labels
  random.treatment <- sample(original.treatment.fem)
  
  #create tabix file
  myobj.all.fem=methRead(file.list.all.fem,
                         sample.id=list("ST2AC10F","ST2AC11F","ST2AC13F","ST2AC14F","ST2AC15F",
                                        "ST2AC16F","ST2AC1F","ST2AC2F","ST2AC3F","ST2AC4F",
                                        "ST2AC5F","ST2AC6F","ST2AC7F","ST2AC8F","ST2AC9F",
                                        "ST2C10F","ST2C11F","ST2C12F","ST2C13F","ST2C14F",
                                        "ST2C15F","ST2C1F","ST2C2F","ST2C3F","ST2C4F",
                                        "ST2C5F","ST2C6F","ST2C7F","ST2C8F","ST2C9F"),
                         assembly="guppyWGBS_shortterm",
                         pipeline="bismarkCoverage",
                         treatment=random.treatment,
                         context="CpG",
                         mincov = 1,
                         dbtype = "tabix",
                         dbdir = "shortterm_allfem_DB"
  )
  
  #run dmr test
  num.DMR.fem[i] <- run.DMR(myobj.all.fem,NULL,keep.chr.allchr,15L)
  
}

# #run permutation test on males 
# for (i in 1:nsim){
#   
#   #randomly shuffle the treatment labels
#   random.treatment <- sample(original.treatment.mal)
#   
#   #create tabix file
#   myobj.all.mal=methRead(file.list.all.mal,
#                          sample.id=list("ST2AC10M","ST2AC11M","ST2AC13M","ST2AC14M","ST2AC15M",
#                                         "ST2AC16M","ST2AC1M","ST2AC2M","ST2AC3M","ST2AC4M",
#                                         "ST2AC5M","ST2AC6M","ST2AC7M","ST2AC8M","ST2AC9M",
#                                         "ST2C10M","ST2C11M","ST2C12M","ST2C13M","ST2C14M",
#                                         "ST2C15M","ST2C1M","ST2C2M","ST2C3M","ST2C4M",
#                                         "ST2C5M","ST2C6M","ST2C7M","ST2C8M", "ST2C9M"),
#                          assembly="guppyWGBS_shortterm",
#                          pipeline="bismarkCoverage",
#                          treatment=random.treatment,
#                          context="CpG",
#                          mincov = 1,
#                          dbtype = "tabix",
#                          dbdir = "shortterm_allmal_DB"
#   )
#   
#   #run dmr test
#   num.DMR.mal[i] <- run.DMR(myobj.all.mal,NULL,keep.chr.allchr,15L)
#   
# }

#save permutation results
saveRDS(num.DMR.fem, file = "./perm_DMR_fem.RDS")
# saveRDS(num.DMR.mal, file = "./perm_DMR_mal.RDS")


