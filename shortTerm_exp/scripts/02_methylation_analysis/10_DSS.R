#####################################################################################################################
### Goal: Convert methylkit obj to BSseq obj and use DSS to call DMS and DMRs
### Author: Janay Fox
### R script
#####################################################################################################################

## Set up ## 

#load in packages
library(DSS)
require(bsseq)
library(methylKit)
library(GenomicRanges)
library(ggplot2)
library(ggpubr)

# read in data 
meth_fem <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit_all/DMSmeth_all_fem_5X.RDS")
meth_mal <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit_all/DMSmeth_all_mal_5X.RDS")
meth_fem_05h <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit/DMSmeth_05h_fem_5X.RDS")
meth_mal_05h <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit/DMSmeth_05h_mal_5X.RDS")
meth_fem_1h <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit/DMSmeth_1h_fem_5X.RDS")
meth_mal_1h <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit/DMSmeth_1h_mal_5X.RDS")
meth_fem_4h <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit/DMSmeth_4h_fem_5X.RDS")
meth_mal_4h <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit/DMSmeth_4h_mal_5X.RDS")
meth_fem_24h <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit/DMSmeth_24h_fem_5X.RDS")
meth_mal_24h <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit/DMSmeth_24h_mal_5X.RDS")
meth_fem_72h <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit/DMSmeth_72h_fem_5X.RDS")
meth_mal_72h <- readRDS("./gup_cue_exp/shortTerm_exp/data/clean/methylKit/DMSmeth_72h_mal_5X.RDS")

#create sample IDs
sampleIDS_fem <- c("ST2AC10F","ST2AC11F","ST2AC13F","ST2AC14F","ST2AC15F",
                   "ST2AC16F","ST2AC1F","ST2AC2F","ST2AC3F","ST2AC4F",
                   "ST2AC5F","ST2AC6F","ST2AC7F","ST2AC8F","ST2AC9F",
                   "ST2C10F","ST2C11F","ST2C12F","ST2C13F","ST2C14F",
                   "ST2C15F","ST2C1F","ST2C2F","ST2C3F","ST2C4F",
                   "ST2C5F","ST2C6F","ST2C7F","ST2C8F","ST2C9F")  

sampleIDS_mal <- c("ST2AC10M","ST2AC11M","ST2AC13M","ST2AC14M","ST2AC15M",
                   "ST2AC16M","ST2AC1M","ST2AC2M","ST2AC3M","ST2AC4M",
                   "ST2AC5M","ST2AC6M","ST2AC7M","ST2AC8M","ST2AC9M",
                   "ST2C10M","ST2C11M","ST2C12M","ST2C13M","ST2C14M",
                   "ST2C15M","ST2C1M","ST2C2M","ST2C3M","ST2C4M",
                   "ST2C5M","ST2C6M","ST2C7M","ST2C8M", "ST2C9M")  

sampleIDS_fem_05h <- c("ST2AC10F","ST2AC16F","ST2AC3F",
                       "ST2C10F","ST2C13F","ST2C2F")

sampleIDS_mal_05h <- c("ST2AC10M","ST2AC16M","ST2AC3M",
                       "ST2C10M","ST2C13M","ST2C2M")

sampleIDS_fem_1h <- c("ST2AC14F","ST2AC2F","ST2AC9F",
                      "ST2C14F","ST2C3F","ST2C9F")

sampleIDS_mal_1h <- c("ST2AC14M","ST2AC2M","ST2AC9M",
                      "ST2C14M","ST2C3M","ST2C9M")

sampleIDS_fem_4h <- c("ST2AC13F","ST2AC1F","ST2AC7F",
                      "ST2C12F","ST2C1F","ST2C7F")

sampleIDS_mal_4h <- c("ST2AC13M","ST2AC1M","ST2AC7M",
                      "ST2C12M","ST2C1M","ST2C7M")

sampleIDS_fem_24h <- c("ST2AC15F","ST2AC4F","ST2AC8F",
                       "ST2C15F","ST2C4F","ST2C8F")

sampleIDS_mal_24h <- c("ST2AC15M","ST2AC4M","ST2AC8M",
                       "ST2C15M","ST2C4M","ST2C8M")

sampleIDS_fem_72h <- c("ST2AC11F","ST2AC5F","ST2AC6F",
                       "ST2C11F","ST2C5F","ST2C6F")

sampleIDS_mal_72h <- c("ST2AC11M","ST2AC5M","ST2AC6M",
                       "ST2C11M","ST2C5M","ST2C6M")

#create groups 
group1_fem = c("ST2AC10F","ST2AC11F","ST2AC13F","ST2AC14F","ST2AC15F",
           "ST2AC16F","ST2AC1F","ST2AC2F","ST2AC3F","ST2AC4F",
           "ST2AC5F","ST2AC6F","ST2AC7F","ST2AC8F","ST2AC9F")
group2_fem = c("ST2C10F","ST2C11F","ST2C12F","ST2C13F","ST2C14F",
           "ST2C15F","ST2C1F","ST2C2F","ST2C3F","ST2C4F",
           "ST2C5F","ST2C6F","ST2C7F","ST2C8F","ST2C9F")

group1_mal = c("ST2AC10M","ST2AC11M","ST2AC13M","ST2AC14M","ST2AC15M",
           "ST2AC16M","ST2AC1M","ST2AC2M","ST2AC3M","ST2AC4M",
           "ST2AC5M","ST2AC6M","ST2AC7M","ST2AC8M","ST2AC9M")
group2_mal = c("ST2C10M","ST2C11M","ST2C12M","ST2C13M","ST2C14M",
           "ST2C15M","ST2C1M","ST2C2M","ST2C3M","ST2C4M",
           "ST2C5M","ST2C6M","ST2C7M","ST2C8M", "ST2C9M")

group1_fem_05h = c("ST2AC10F","ST2AC16F","ST2AC3F")
group2_fem_05h = c("ST2C10F","ST2C13F","ST2C2F")

group1_mal_05h = c("ST2AC10M","ST2AC16M","ST2AC3M")
group2_mal_05h = c("ST2C10M","ST2C13M","ST2C2M")

group1_fem_1h = c("ST2AC14F","ST2AC2F","ST2AC9F")
group2_fem_1h = c("ST2C14F","ST2C3F","ST2C9F")

group1_mal_1h = c("ST2AC14M","ST2AC2M","ST2AC9M")
group2_mal_1h = c("ST2C14M","ST2C3M","ST2C9M")

group1_fem_4h = c("ST2AC13F","ST2AC1F","ST2AC7F")
group2_fem_4h = c("ST2C12F","ST2C1F","ST2C7F")

group1_mal_4h = c("ST2AC13M","ST2AC1M","ST2AC7M")
group2_mal_4h = c("ST2C12M","ST2C1M","ST2C7M")

group1_fem_24h = c("ST2AC15F","ST2AC4F","ST2AC8F")
group2_fem_24h = c("ST2C15F","ST2C4F","ST2C8F")

group1_mal_24h = c("ST2AC15M","ST2AC4M","ST2AC8M")
group2_mal_24h = c("ST2C15M","ST2C4M","ST2C8M")

group1_fem_72h = c("ST2AC11F","ST2AC5F","ST2AC6F")
group2_fem_72h = c("ST2C11F","ST2C5F","ST2C6F")

group1_mal_72h = c("ST2AC11M","ST2AC5M","ST2AC6M")
group2_mal_72h = c("ST2C11M","ST2C5M","ST2C6M")

## Create functions ##
#create function to run DSS on each time point 
methylKitToDSStp <- function(methobj, sampleIDs, Group1, Group2){
  
  #create BSseq obj
  rowRanges <- GRanges(
    seqnames = methobj$chr,
    ranges   = IRanges(start = methobj$start, end = methobj$end),
    strand   = methobj$strand
  )
  
  
  Cov <- cbind(methobj$coverage1, methobj$coverage2, methobj$coverage3,
               methobj$coverage4, methobj$coverage5, methobj$coverage6
  )
  
  M <- cbind(methobj$numCs1, methobj$numCs2, methobj$numCs3,
             methobj$numCs4, methobj$numCs5, methobj$numCs6
  )
  
  colnames(Cov) <- sampleIDs
  colnames(M)   <- sampleIDs
  
  BSobj <- BSseq(
    M = M,
    Cov = Cov,
    gr = rowRanges
  )
  
  #do DMS tests 
  dmlTest = DMLtest(BSobj, group1 = Group1,
                        group2 = Group2,
                        smoothing=TRUE)
  
  #call DMSs
  dmcs = callDML(dmlTest, delta = 0.2, p.threshold=0.0125)
  
  #call DMRs
  dmrs <- callDMR(dmlTest, p.threshold = 0.0125, delta = 0.2, minlen = 100, minCG = 5)
  
  results <- list(dmlTest = dmlTest,
                  dmcs = dmcs,
                  dmrs = dmrs)
  
  return(results)
}

#create function to run DSS on all time points 
methylKitToDSSall <- function(methobj, sampleIDs, Group1, Group2){
  
  #create BSseq obj
  rowRanges <- GRanges(
    seqnames = methobj$chr,
    ranges   = IRanges(start = methobj$start, end = methobj$end),
    strand   = methobj$strand
  )
  
  
  Cov <- cbind(methobj$coverage1, methobj$coverage2, methobj$coverage3,
               methobj$coverage4, methobj$coverage5, methobj$coverage6,
               methobj$coverage7, methobj$coverage8, methobj$coverage9,
               methobj$coverage10, methobj$coverage11, methobj$coverage12,
               methobj$coverage13, methobj$coverage14, methobj$coverage15,
               methobj$coverage16, methobj$coverage17, methobj$coverage18,
               methobj$coverage19, methobj$coverage20, methobj$coverage21,
               methobj$coverage22, methobj$coverage23, methobj$coverage24,
               methobj$coverage25, methobj$coverage26, methobj$coverage27,
               methobj$coverage28, methobj$coverage29, methobj$coverage30
  )
  
  M <- cbind(methobj$numCs1, methobj$numCs2, methobj$numCs3,
             methobj$numCs4, methobj$numCs5, methobj$numCs6,
             methobj$numCs7, methobj$numCs8, methobj$numCs9,
             methobj$numCs10, methobj$numCs11, methobj$numCs12,
             methobj$numCs13, methobj$numCs14, methobj$numCs15,
             methobj$numCs16, methobj$numCs17, methobj$numCs18,
             methobj$numCs19, methobj$numCs20, methobj$numCs21,
             methobj$numCs22, methobj$numCs23, methobj$numCs24,
             methobj$numCs25, methobj$numCs26, methobj$numCs27,
             methobj$numCs28, methobj$numCs29, methobj$numCs30
  )
  
  colnames(Cov) <- sampleIDs
  colnames(M)   <- sampleIDs
  
  BSobj <- BSseq(
    M = M,
    Cov = Cov,
    gr = rowRanges
  )
  
  #do DMS tests 
  dmlTest = DMLtest(BSobj, group1 = Group1,
                    group2 = Group2,
                    smoothing=TRUE)
  
  #call DMSs
  dmcs = callDML(dmlTest, delta = 0.2, p.threshold=0.0125)
  
  #call DMRs
  dmrs <- callDMR(dmlTest, p.threshold = 0.0125, delta = 0.2, minlen = 100, minCG = 5)
  
  results <- list(dmlTest = dmlTest,
                  dmcs = dmcs,
                  dmrs = dmrs)
  
  return(results)
}

#plot barplots 
barplot.funct <- function(bar.data, site, ylabel){
  ggplot(data = bar.data, aes(x=time_point, y = {{site}}, fill = sex)) + 
    geom_bar(stat="identity", position = position_dodge()) +
    theme_bw() + xlab("Time Point") + ylab(ylabel) +
    scale_fill_manual(values = c("fem" = "#481567FF", "mal" = "#9AC8CD"), name = "Sex", 
                      labels = c("fem" = "Female", "mal" = "Male")) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text = element_text(size = 11),
          legend.title = element_text(size = 12)) 
}

## Run DSS ## 
DSS_fem_05h <- methylKitToDSStp(meth_fem_05h, sampleIDS_fem_05h, group1_fem_05h, group2_fem_05h)
DSS_mal_05h <- methylKitToDSStp(meth_mal_05h, sampleIDS_mal_05h, group1_mal_05h, group2_mal_05h)

DSS_fem_1h <- methylKitToDSStp(meth_fem_1h, sampleIDS_fem_1h, group1_fem_1h, group2_fem_1h)
DSS_mal_1h <- methylKitToDSStp(meth_mal_1h, sampleIDS_mal_1h, group1_mal_1h, group2_mal_1h)

DSS_fem_4h <- methylKitToDSStp(meth_fem_4h, sampleIDS_fem_4h, group1_fem_4h, group2_fem_4h)
DSS_mal_4h <- methylKitToDSStp(meth_mal_4h, sampleIDS_mal_4h, group1_mal_4h, group2_mal_4h)

DSS_fem_24h <- methylKitToDSStp(meth_fem_24h, sampleIDS_fem_24h, group1_fem_24h, group2_fem_24h)
DSS_mal_24h <- methylKitToDSStp(meth_mal_24h, sampleIDS_mal_24h, group1_mal_24h, group2_mal_24h)

DSS_fem_72h <- methylKitToDSStp(meth_fem_72h, sampleIDS_fem_72h, group1_fem_72h, group2_fem_72h)
DSS_mal_72h <- methylKitToDSStp(meth_mal_72h, sampleIDS_mal_72h, group1_mal_72h, group2_mal_72h)

DSS_fem_all <- methylKitToDSSall(meth_fem, sampleIDS_fem, group1_fem, group2_fem)
DSS_mal_all <- methylKitToDSSall(meth_mal, sampleIDS_mal, group1_mal, group2_mal)

## Plot ## 

DSS_fem_1h$dmcs <- callDML(DSS_fem_1h$dmlTest, delta = 0.2,p.threshold = 0.0125)

#get number of DMSs and DMRs
DMS.list <- c(nrow(DSS_fem_all$dmcs), nrow(DSS_fem_05h$dmcs), nrow(DSS_fem_1h$dmcs), nrow(DSS_fem_4h$dmcs), 
                 nrow(DSS_fem_24h$dmcs), nrow(DSS_fem_72h$dmcs),
                 nrow(DSS_mal_all$dmcs), nrow(DSS_mal_05h$dmcs), nrow(DSS_mal_1h$dmcs), nrow(DSS_mal_4h$dmcs),
                 nrow(DSS_mal_24h$dmcs), nrow(DSS_mal_72h$dmcs))

DMR.list <- c(0, nrow(DSS_fem_05h$dmrs), nrow(DSS_fem_1h$dmrs), nrow(DSS_fem_4h$dmrs), nrow(DSS_fem_24h$dmrs), 
                 nrow(DSS_fem_72h$dmrs),
                 0, nrow(DSS_mal_05h$dmrs), nrow(DSS_mal_1h$dmrs), nrow(DSS_mal_4h$dmrs), nrow(DSS_mal_24h$dmrs), 
                 nrow(DSS_mal_72h$dmrs))

barplot.data <- data.frame(num_DMS = DMS.list, num_DMR = DMR.list, 
                           time_point = c("all", "0.5h", "1h", "4h", "24h", "72h",
                                          "all", "0.5h", "1h", "4h", "24h", "72h"),
                           sex = c("fem", "fem", "fem", "fem", "fem", "fem",
                                   "mal", "mal", "mal", "mal", "mal", "mal"))

barplot.data$time_point <- factor(barplot.data$time_point, ordered = TRUE, 
                                  levels = c("0.5h", "1h", "4h", "24h", "72h", "all"))

#plot barplots
DMS.barplot <- barplot.funct(barplot.data, num_DMS, "Number of DMSs")
DMS.barplot

DMR.barplot <- barplot.funct(barplot.data, num_DMR, "Number of DMRs")
DMR.barplot

#make a panel 
barplot.panel <- ggarrange(DMS.barplot, DMR.barplot, labels = c("A", "B"),
                           nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom") 
#save plots
tiff("./gup_cue_exp/shortTerm_exp/DSS_panel_barPlot.tiff", units="in", width = 6, height = 4, res = 600)
barplot.panel
dev.off()
