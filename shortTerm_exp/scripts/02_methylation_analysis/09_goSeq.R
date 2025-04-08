#######################################################
### Goal: Run and plot GO-seq analysis
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")

#load packaes
library(ggplot2)
library(dplyr)
library(genomation)
library(DescTools)
library(tibble)
library(GenomicRanges)
library(methylKit)
library(GOstats)
library(biomaRt)
library(GSEABase)
library(ggpubr)
library(stringr)

#write functions 
#function for running goseq
run.gostat <- function(test.data, uni.data, ont){
  #get background gene sets from all CpGs
  uni.data <- getAssociationWithTSS(uni.data)
  
  #limit to genes within 10kb 
  test.data <- subset(test.data, dist.to.feature >= -10000 & dist.to.feature <= 10000)
  
  universe <- uni.data$feature.name
  genes <- test.data$feature.name
  params <- GSEAGOHyperGParams(name = "shorterm GSEA based annot Params",
                               geneSetCollection = gsc,
                               geneIds = genes,
                               universeGeneIds = universe,
                               ontology = ont,
                               pvalueCutoff = 0.05,
                               conditional = FALSE,
                               testDirection = "over")
  
  over <- hyperGTest(params)
  print(over)
  sum.over <- summary(over)
  sum.over$FDR <-  p.adjust(sum.over$Pvalue, method = "fdr")
  return(sum.over)
}

#function for plotting go seq
plot.funct <- function(data.05h, data.1h, data.4h, data.24h, data.72h,
                       total.05h, total.1h, total.4h, total.24h, total.72h,
                       name.plot){
  #calculate percent
  data.05h$perc <- data.05h$Count / total.05h * 100
  data.1h$perc <- data.1h$Count / total.1h * 100
  data.4h$perc <- data.4h$Count / total.4h * 100
  data.24h$perc <- data.24h$Count / total.24h * 100
  data.72h$perc <- data.72h$Count / total.72h * 100
  
  #select top 10 by p value and percent
  data.05h <- data.05h %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  data.1h <- data.1h %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  data.4h <- data.4h %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  data.24h <- data.24h %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  data.72h <- data.72h %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  
  #add time point variable
  data.05h$tp <- "0.5h"
  data.1h$tp <- "1h"
  data.4h$tp <- "4h"
  data.24h$tp <- "24h"
  data.72h$tp <- "72h"
  
  #bind all data
  data <- rbind(data.05h, data.1h, data.4h, data.24h, data.72h)
  
  #convert time point to ordered factor
  data$tp <- factor(data$tp, levels = c("0.5h", "1h", "4h", "24h", "72h"), ordered = TRUE)
  
  #plot
  p <- data %>% ggplot(aes(x = tp, y = Term, size = perc, color = tp, alpha = FDR)) + geom_point() + theme_bw() +
    scale_color_manual(values = c("#B8DE29FF", "#55C667FF", "#1F968BFF", "#39568CFF", "#481567FF")) +
    scale_size(range = c(3,10), name = "% of total", breaks = c(5,10,25), limits = c(0,50)) + ylab(NULL) + xlab("Time Point") +
    scale_alpha_continuous(name = "Adj. P value", limits = c(0,0.05), range = c(0.8,0.1),
                           breaks = c(0.0001, 0.001, 0.01, 0.05), labels = c("0.0001", 0.001, 0.01, 0.05)) + guides(color="none") +
    scale_y_discrete(labels = function(Term) str_wrap(Term, width = 100)) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 10), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
  
  ggsave(filename = name.plot, plot = p, width = 8, height = 10, units = "in", dpi = 300)
  
  return(p)
}

#function for plotting top 5
plot.funct.top5 <- function(data.05h, data.1h, data.4h, data.24h, data.72h,
                            total.05h, total.1h, total.4h, total.24h, total.72h,
                            name.plot){
  #calculate percent
  data.05h$perc <- data.05h$Count / total.05h * 100
  data.1h$perc <- data.1h$Count / total.1h * 100
  data.4h$perc <- data.4h$Count / total.4h * 100
  data.24h$perc <- data.24h$Count / total.24h * 100
  data.72h$perc <- data.72h$Count / total.72h * 100
  
  #select top 10 by p value and percent
  data.05h <- data.05h %>% slice_min(tibble(FDR, perc), n = 5, with_ties = FALSE)
  data.1h <- data.1h %>% slice_min(tibble(FDR, perc), n = 5, with_ties = FALSE)
  data.4h <- data.4h %>% slice_min(tibble(FDR, perc), n = 5, with_ties = FALSE)
  data.24h <- data.24h %>% slice_min(tibble(FDR, perc), n = 5, with_ties = FALSE)
  data.72h <- data.72h %>% slice_min(tibble(FDR, perc), n = 5, with_ties = FALSE)
  
  #add time point variable
  data.05h$tp <- "0.5h"
  data.1h$tp <- "1h"
  data.4h$tp <- "4h"
  data.24h$tp <- "24h"
  data.72h$tp <- "72h"
  
  #bind all data
  data <- rbind(data.05h, data.1h, data.4h, data.24h, data.72h)
  
  #convert time point to ordered factor
  data$tp <- factor(data$tp, levels = c("0.5h", "1h", "4h", "24h", "72h"), ordered = TRUE)
  
  #plot
  p <- data %>% ggplot(aes(x = tp, y = Term, size = perc, color = tp, alpha = FDR)) + geom_point() + theme_bw() +
    scale_color_manual(values = c("#B8DE29FF", "#55C667FF", "#1F968BFF", "#39568CFF", "#481567FF")) +
    scale_size(range = c(3,10), name = "% of total", breaks = c(5,10,25), limits = c(0,50)) + ylab(NULL) + xlab("Time Point") +
    scale_alpha_continuous(name = "Adj. P value", limits = c(0,0.05), range = c(0.8,0.1),
                           breaks = c(0.0001, 0.001, 0.01, 0.05), labels = c("0.0001", 0.001, 0.01, 0.05)) + guides(color="none") +
    scale_y_discrete(labels = function(Term) str_wrap(Term, width = 60)) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 10), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
  
  ggsave(filename = name.plot, plot = p, width = 9, height = 10, units = "in", dpi = 300)
  
  return(p)
}

#function for getting genes
gene.list.funct <- function(data, name.csv){
  #limit to genes within 10kb
  data <- subset(data, dist.to.feature >= -10000 & dist.to.feature <= 10000)
  #get list to seach in biomart
  val <- data$feature.name
  #search biomart
  gene.data <- getBM(mart = bm, attributes = c('ensembl_transcript_id','external_gene_name','go_id', "name_1006",
                                               "namespace_1003", "go_linkage_type"), filters = "ensembl_transcript_id", values = val)
  gene.data <- gene.data[!duplicated(gene.data),]
  write.csv(gene.data, name.csv, row.names = FALSE)
}

#function for plotting pooled analysis 
plot_pooled <- function(data.hyper, total.hyper, data.hypo, total.hypo, name.plot){
  #add labels
  data.hyper$direction <- "hyper"
  data.hypo$direction <- "hypo"
  
  #calculate percent
  data.hyper$perc <- data.hyper$Count / total.hyper * 100
  data.hypo$perc <- data.hypo$Count / total.hypo * 100
  
  #select top 10 by p value and percent
  data.hyper <- data.hyper %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  data.hypo <- data.hypo %>% slice_min(tibble(FDR, perc), n = 10, with_ties = FALSE)
  
  #merge data
  data <- rbind(data.hyper, data.hypo)
  
  #plot
  p <- data %>% ggplot(aes(x = direction, y = Term, color = direction, size = perc, alpha = FDR)) + geom_point() + theme_bw() +
    scale_color_manual(values = c("#55C667FF", "#39568CFF")) +
    scale_size(range = c(3,10), name = "% of total", breaks = c(5,10,25), limits = c(0,50)) + ylab(NULL) + xlab("Time Point") +
    scale_alpha_continuous(name = "Adj. P value", limits = c(0,0.05), range = c(0.8,0.1),
                           breaks = c(0.0001, 0.001, 0.01, 0.05), labels = c("0.0001", 0.001, 0.01, 0.05)) + guides(color="none") +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 10), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
  
  ggsave(filename = name.plot, plot = p, width = 8, height = 10, units = "in", dpi = 300)
}

#load in data
DMS_TSS_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMS_TSS_all_fem.RDS")
DMS_TSS_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMS_TSS_05h_fem.RDS")
DMS_TSS_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMS_TSS_1h_fem.RDS")
DMS_TSS_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMS_TSS_4h_fem.RDS")
DMS_TSS_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMS_TSS_24h_fem.RDS")
DMS_TSS_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMS_TSS_72h_fem.RDS")

DMS_TSS_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMS_TSS_all_mal.RDS")
DMS_TSS_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMS_TSS_05h_mal.RDS")
DMS_TSS_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMS_TSS_1h_mal.RDS")
DMS_TSS_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMS_TSS_4h_mal.RDS")
DMS_TSS_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMS_TSS_24h_mal.RDS")
DMS_TSS_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMS_TSS_72h_mal.RDS")

DMR_TSS_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMR_TSS_all_fem.RDS")
DMR_TSS_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMR_TSS_05h_fem.RDS")
DMR_TSS_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMR_TSS_1h_fem.RDS")
DMR_TSS_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMR_TSS_4h_fem.RDS")
DMR_TSS_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMR_TSS_24h_fem.RDS")
DMR_TSS_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMR_TSS_72h_fem.RDS")

DMR_TSS_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMR_TSS_all_mal.RDS")
DMR_TSS_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMR_TSS_05h_mal.RDS")
DMR_TSS_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMR_TSS_1h_mal.RDS")
DMR_TSS_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMR_TSS_4h_mal.RDS")
DMR_TSS_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMR_TSS_24h_mal.RDS")
DMR_TSS_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMR_TSS_72h_mal.RDS")

DMS_TSS_all_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMS_TSS_all_fem_hyper.RDS")
DMS_TSS_05h_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMS_TSS_05h_fem_hyper.RDS")
DMS_TSS_1h_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMS_TSS_1h_fem_hyper.RDS")
DMS_TSS_4h_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMS_TSS_4h_fem_hyper.RDS")
DMS_TSS_24h_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMS_TSS_24h_fem_hyper.RDS")
DMS_TSS_72h_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMS_TSS_72h_fem_hyper.RDS")

DMS_TSS_all_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMS_TSS_all_mal_hyper.RDS")
DMS_TSS_05h_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMS_TSS_05h_mal_hyper.RDS")
DMS_TSS_1h_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMS_TSS_1h_mal_hyper.RDS")
DMS_TSS_4h_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMS_TSS_4h_mal_hyper.RDS")
DMS_TSS_24h_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMS_TSS_24h_mal_hyper.RDS")
DMS_TSS_72h_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMS_TSS_72h_mal_hyper.RDS")

DMR_TSS_all_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMR_TSS_all_fem_hyper.RDS")
DMR_TSS_05h_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMR_TSS_05h_fem_hyper.RDS")
DMR_TSS_1h_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMR_TSS_1h_fem_hyper.RDS")
DMR_TSS_4h_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMR_TSS_4h_fem_hyper.RDS")
DMR_TSS_24h_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMR_TSS_24h_fem_hyper.RDS")
DMR_TSS_72h_fem_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMR_TSS_72h_fem_hyper.RDS")

DMR_TSS_all_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMR_TSS_all_mal_hyper.RDS")
DMR_TSS_05h_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMR_TSS_05h_mal_hyper.RDS")
DMR_TSS_1h_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMR_TSS_1h_mal_hyper.RDS")
DMR_TSS_4h_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMR_TSS_4h_mal_hyper.RDS")
DMR_TSS_24h_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMR_TSS_24h_mal_hyper.RDS")
DMR_TSS_72h_mal_hyper <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMR_TSS_72h_mal_hyper.RDS")

DMS_TSS_all_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMS_TSS_all_fem_hypo.RDS")
DMS_TSS_05h_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMS_TSS_05h_fem_hypo.RDS")
DMS_TSS_1h_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMS_TSS_1h_fem_hypo.RDS")
DMS_TSS_4h_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMS_TSS_4h_fem_hypo.RDS")
DMS_TSS_24h_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMS_TSS_24h_fem_hypo.RDS")
DMS_TSS_72h_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMS_TSS_72h_fem_hypo.RDS")

DMS_TSS_all_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMS_TSS_all_mal_hypo.RDS")
DMS_TSS_05h_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMS_TSS_05h_mal_hypo.RDS")
DMS_TSS_1h_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMS_TSS_1h_mal_hypo.RDS")
DMS_TSS_4h_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMS_TSS_4h_mal_hypo.RDS")
DMS_TSS_24h_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMS_TSS_24h_mal_hypo.RDS")
DMS_TSS_72h_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMS_TSS_72h_mal_hypo.RDS")

DMR_TSS_all_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMR_TSS_all_fem_hypo.RDS")
DMR_TSS_05h_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMR_TSS_05h_fem_hypo.RDS")
DMR_TSS_1h_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMR_TSS_1h_fem_hypo.RDS")
DMR_TSS_4h_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMR_TSS_4h_fem_hypo.RDS")
DMR_TSS_24h_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMR_TSS_24h_fem_hypo.RDS")
DMR_TSS_72h_fem_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMR_TSS_72h_fem_hypo.RDS")

DMR_TSS_all_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/DMR_TSS_all_mal_hypo.RDS")
DMR_TSS_05h_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/DMR_TSS_05h_mal_hypo.RDS")
DMR_TSS_1h_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/DMR_TSS_1h_mal_hypo.RDS")
DMR_TSS_4h_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/DMR_TSS_4h_mal_hypo.RDS")
DMR_TSS_24h_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/DMR_TSS_24h_mal_hypo.RDS")
DMR_TSS_72h_mal_hypo <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/DMR_TSS_72h_mal_hypo.RDS")

CpG_anno_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/CpG_anno_all_fem.RDS")
CpG_anno_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/CpG_anno_05h_fem.RDS")
CpG_anno_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/CpG_anno_1h_fem.RDS")
CpG_anno_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/CpG_anno_4h_fem.RDS")
CpG_anno_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/CpG_anno_24h_fem.RDS")
CpG_anno_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/CpG_anno_72h_fem.RDS")

CpG_anno_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/CpG_anno_all_mal.RDS")
CpG_anno_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/CpG_anno_05h_mal.RDS")
CpG_anno_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/CpG_anno_1h_mal.RDS")
CpG_anno_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/CpG_anno_4h_mal.RDS")
CpG_anno_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/CpG_anno_24h_mal.RDS")
CpG_anno_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/CpG_anno_72h_mal.RDS")

regions_anno_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/regions_anno_all_fem.RDS")
regions_anno_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/regions_anno_05h_fem.RDS")
regions_anno_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/regions_anno_1h_fem.RDS")
regions_anno_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/regions_anno_4h_fem.RDS")
regions_anno_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/regions_anno_24h_fem.RDS")
regions_anno_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/regions_anno_72h_fem.RDS")

regions_anno_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_all/regions_anno_all_mal.RDS")
regions_anno_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_05h/regions_anno_05h_mal.RDS")
regions_anno_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_1h/regions_anno_1h_mal.RDS")
regions_anno_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_4h/regions_anno_4h_mal.RDS")
regions_anno_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_24h/regions_anno_24h_mal.RDS")
regions_anno_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/anno_res_72h/regions_anno_72h_mal.RDS")

## Go seq analysis ##
#get proper biomart
bm <- useEnsembl(biomart = "ensembl")
bm <- useDataset("preticulata_gene_ensembl", mart = bm)

#grab go terms
go_list <- getBM(mart = bm, attributes = c('ensembl_transcript_id','external_gene_name','go_id', "name_1006", 
                                                                              "namespace_1003", "go_linkage_type"))
#prep go terms for gene mappings
go_list_data <- go_list[,c(3,6,1)] #limit to needed data
go_list_data <- go_list_data[go_list_data$go_id != '',] #remove empty go terms

#set go frame
goFrame <- GOFrame(go_list_data, organism = "guppy")
goAllFrame <- GOAllFrame(goFrame)
#set gene set
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

#run go seg on all gene sets
DMS_all_fem_gostat_BP <- run.gostat(DMS_TSS_all_fem, CpG_anno_all_fem, "BP")
DMS_05h_fem_gostat_BP <- run.gostat(DMS_TSS_05h_fem, CpG_anno_05h_fem, "BP")
DMS_1h_fem_gostat_BP <- run.gostat(DMS_TSS_1h_fem, CpG_anno_1h_fem, "BP")
DMS_4h_fem_gostat_BP<- run.gostat(DMS_TSS_4h_fem, CpG_anno_4h_fem, "BP")
DMS_24h_fem_gostat_BP <- run.gostat(DMS_TSS_24h_fem, CpG_anno_24h_fem, "BP")
DMS_72h_fem_gostat_BP <- run.gostat(DMS_TSS_72h_fem, CpG_anno_72h_fem, "BP")

DMS_all_mal_gostat_BP <- run.gostat(DMS_TSS_all_mal, CpG_anno_all_mal, "BP")
DMS_05h_mal_gostat_BP <- run.gostat(DMS_TSS_05h_mal, CpG_anno_05h_mal, "BP")
DMS_1h_mal_gostat_BP <- run.gostat(DMS_TSS_1h_mal, CpG_anno_1h_mal, "BP")
DMS_4h_mal_gostat_BP <- run.gostat(DMS_TSS_4h_mal, CpG_anno_4h_mal, "BP")
DMS_24h_mal_gostat_BP <- run.gostat(DMS_TSS_24h_mal, CpG_anno_24h_mal, "BP")
DMS_72h_mal_gostat_BP <- run.gostat(DMS_TSS_72h_mal, CpG_anno_72h_mal, "BP")

DMS_all_fem_gostat_MF <- run.gostat(DMS_TSS_all_fem, CpG_anno_all_fem, "MF")
DMS_05h_fem_gostat_MF <- run.gostat(DMS_TSS_05h_fem, CpG_anno_05h_fem, "MF")
DMS_1h_fem_gostat_MF <- run.gostat(DMS_TSS_1h_fem, CpG_anno_1h_fem, "MF")
DMS_4h_fem_gostat_MF<- run.gostat(DMS_TSS_4h_fem, CpG_anno_4h_fem, "MF")
DMS_24h_fem_gostat_MF <- run.gostat(DMS_TSS_24h_fem, CpG_anno_24h_fem, "MF")
DMS_72h_fem_gostat_MF <- run.gostat(DMS_TSS_72h_fem, CpG_anno_72h_fem, "MF")

DMS_all_mal_gostat_MF <- run.gostat(DMS_TSS_all_mal, CpG_anno_all_mal, "MF")
DMS_05h_mal_gostat_MF <- run.gostat(DMS_TSS_05h_mal, CpG_anno_05h_mal, "MF")
DMS_1h_mal_gostat_MF <- run.gostat(DMS_TSS_1h_mal, CpG_anno_1h_mal, "MF")
DMS_4h_mal_gostat_MF <- run.gostat(DMS_TSS_4h_mal, CpG_anno_4h_mal, "MF")
DMS_24h_mal_gostat_MF <- run.gostat(DMS_TSS_24h_mal, CpG_anno_24h_mal, "MF")
DMS_72h_mal_gostat_MF <- run.gostat(DMS_TSS_72h_mal, CpG_anno_72h_mal, "MF")

DMS_all_fem_gostat_CC <- run.gostat(DMS_TSS_all_fem, CpG_anno_all_fem, "CC")
DMS_05h_fem_gostat_CC <- run.gostat(DMS_TSS_05h_fem, CpG_anno_05h_fem, "CC")
DMS_1h_fem_gostat_CC <- run.gostat(DMS_TSS_1h_fem, CpG_anno_1h_fem, "CC")
DMS_4h_fem_gostat_CC<- run.gostat(DMS_TSS_4h_fem, CpG_anno_4h_fem, "CC")
DMS_24h_fem_gostat_CC <- run.gostat(DMS_TSS_24h_fem, CpG_anno_24h_fem, "CC")
DMS_72h_fem_gostat_CC <- run.gostat(DMS_TSS_72h_fem, CpG_anno_72h_fem, "CC")

DMS_all_mal_gostat_CC <- run.gostat(DMS_TSS_all_mal, CpG_anno_all_mal, "CC")
DMS_05h_mal_gostat_CC <- run.gostat(DMS_TSS_05h_mal, CpG_anno_05h_mal, "CC")
DMS_1h_mal_gostat_CC <- run.gostat(DMS_TSS_1h_mal, CpG_anno_1h_mal, "CC")
DMS_4h_mal_gostat_CC <- run.gostat(DMS_TSS_4h_mal, CpG_anno_4h_mal, "CC")
DMS_24h_mal_gostat_CC <- run.gostat(DMS_TSS_24h_mal, CpG_anno_24h_mal, "CC")
DMS_72h_mal_gostat_CC <- run.gostat(DMS_TSS_72h_mal, CpG_anno_72h_mal, "CC")

DMR_all_fem_gostat_BP <- run.gostat(DMR_TSS_all_fem, regions_anno_all_fem, "BP")
DMR_05h_fem_gostat_BP <- run.gostat(DMR_TSS_05h_fem, regions_anno_05h_fem, "BP")
DMR_1h_fem_gostat_BP <- run.gostat(DMR_TSS_1h_fem, regions_anno_1h_fem, "BP")
DMR_4h_fem_gostat_BP<- run.gostat(DMR_TSS_4h_fem, regions_anno_4h_fem, "BP")
DMR_24h_fem_gostat_BP <- run.gostat(DMR_TSS_24h_fem, regions_anno_24h_fem, "BP")
DMR_72h_fem_gostat_BP <- run.gostat(DMR_TSS_72h_fem, regions_anno_72h_fem, "BP")

DMR_all_mal_gostat_BP <- run.gostat(DMR_TSS_all_mal, regions_anno_all_mal, "BP")
DMR_05h_mal_gostat_BP <- run.gostat(DMR_TSS_05h_mal, regions_anno_05h_mal, "BP")
DMR_1h_mal_gostat_BP <- run.gostat(DMR_TSS_1h_mal, regions_anno_1h_mal, "BP")
DMR_4h_mal_gostat_BP <- run.gostat(DMR_TSS_4h_mal, regions_anno_4h_mal, "BP")
DMR_24h_mal_gostat_BP <- run.gostat(DMR_TSS_24h_mal, regions_anno_24h_mal, "BP")
DMR_72h_mal_gostat_BP <- run.gostat(DMR_TSS_72h_mal, regions_anno_72h_mal, "BP")

DMR_all_fem_gostat_MF <- run.gostat(DMR_TSS_all_fem, regions_anno_all_fem, "MF")
DMR_05h_fem_gostat_MF <- run.gostat(DMR_TSS_05h_fem, regions_anno_05h_fem, "MF")
DMR_1h_fem_gostat_MF <- run.gostat(DMR_TSS_1h_fem, regions_anno_1h_fem, "MF")
DMR_4h_fem_gostat_MF<- run.gostat(DMR_TSS_4h_fem, regions_anno_4h_fem, "MF")
DMR_24h_fem_gostat_MF <- run.gostat(DMR_TSS_24h_fem, regions_anno_24h_fem, "MF")
DMR_72h_fem_gostat_MF <- run.gostat(DMR_TSS_72h_fem, regions_anno_72h_fem, "MF")

DMR_all_mal_gostat_MF <- run.gostat(DMR_TSS_all_mal, regions_anno_all_mal, "MF")
DMR_05h_mal_gostat_MF <- run.gostat(DMR_TSS_05h_mal, regions_anno_05h_mal, "MF")
DMR_1h_mal_gostat_MF <- run.gostat(DMR_TSS_1h_mal, regions_anno_1h_mal, "MF")
DMR_4h_mal_gostat_MF <- run.gostat(DMR_TSS_4h_mal, regions_anno_4h_mal, "MF")
DMR_24h_mal_gostat_MF <- run.gostat(DMR_TSS_24h_mal, regions_anno_24h_mal, "MF")
DMR_72h_mal_gostat_MF <- run.gostat(DMR_TSS_72h_mal, regions_anno_72h_mal, "MF")

DMR_all_fem_gostat_CC <- run.gostat(DMR_TSS_all_fem, regions_anno_all_fem, "CC")
DMR_05h_fem_gostat_CC <- run.gostat(DMR_TSS_05h_fem, regions_anno_05h_fem, "CC")
DMR_1h_fem_gostat_CC <- run.gostat(DMR_TSS_1h_fem, regions_anno_1h_fem, "CC")
DMR_4h_fem_gostat_CC<- run.gostat(DMR_TSS_4h_fem, regions_anno_4h_fem, "CC")
DMR_24h_fem_gostat_CC <- run.gostat(DMR_TSS_24h_fem, regions_anno_24h_fem, "CC")
DMR_72h_fem_gostat_CC <- run.gostat(DMR_TSS_72h_fem, regions_anno_72h_fem, "CC")

DMR_all_mal_gostat_CC <- run.gostat(DMR_TSS_all_mal, regions_anno_all_mal, "CC")
DMR_05h_mal_gostat_CC <- run.gostat(DMR_TSS_05h_mal, regions_anno_05h_mal, "CC")
DMR_1h_mal_gostat_CC <- run.gostat(DMR_TSS_1h_mal, regions_anno_1h_mal, "CC")
DMR_4h_mal_gostat_CC <- run.gostat(DMR_TSS_4h_mal, regions_anno_4h_mal, "CC")
DMR_24h_mal_gostat_CC <- run.gostat(DMR_TSS_24h_mal, regions_anno_24h_mal, "CC")
DMR_72h_mal_gostat_CC <- run.gostat(DMR_TSS_72h_mal, regions_anno_72h_mal, "CC")

#run on hyper gene sets 
DMS_all_fem_gostat_BP_hyper <- run.gostat(DMS_TSS_all_fem_hyper, CpG_anno_all_fem, "BP")
DMS_05h_fem_gostat_BP_hyper <- run.gostat(DMS_TSS_05h_fem_hyper, CpG_anno_05h_fem, "BP")
DMS_1h_fem_gostat_BP_hyper <- run.gostat(DMS_TSS_1h_fem_hyper, CpG_anno_1h_fem, "BP")
DMS_4h_fem_gostat_BP_hyper<- run.gostat(DMS_TSS_4h_fem_hyper, CpG_anno_4h_fem, "BP")
DMS_24h_fem_gostat_BP_hyper <- run.gostat(DMS_TSS_24h_fem_hyper, CpG_anno_24h_fem, "BP")
DMS_72h_fem_gostat_BP_hyper <- run.gostat(DMS_TSS_72h_fem_hyper, CpG_anno_72h_fem, "BP")

DMS_all_mal_gostat_BP_hyper <- run.gostat(DMS_TSS_all_mal_hyper, CpG_anno_all_mal, "BP")
DMS_05h_mal_gostat_BP_hyper <- run.gostat(DMS_TSS_05h_mal_hyper, CpG_anno_05h_mal, "BP")
DMS_1h_mal_gostat_BP_hyper <- run.gostat(DMS_TSS_1h_mal_hyper, CpG_anno_1h_mal, "BP")
DMS_4h_mal_gostat_BP_hyper <- run.gostat(DMS_TSS_4h_mal_hyper, CpG_anno_4h_mal, "BP")
DMS_24h_mal_gostat_BP_hyper <- run.gostat(DMS_TSS_24h_mal_hyper, CpG_anno_24h_mal, "BP")
DMS_72h_mal_gostat_BP_hyper <- run.gostat(DMS_TSS_72h_mal_hyper, CpG_anno_72h_mal, "BP")

DMS_all_fem_gostat_MF_hyper <- run.gostat(DMS_TSS_all_fem_hyper, CpG_anno_all_fem, "MF")
DMS_05h_fem_gostat_MF_hyper <- run.gostat(DMS_TSS_05h_fem_hyper, CpG_anno_05h_fem, "MF")
DMS_1h_fem_gostat_MF_hyper <- run.gostat(DMS_TSS_1h_fem_hyper, CpG_anno_1h_fem, "MF")
DMS_4h_fem_gostat_MF_hyper<- run.gostat(DMS_TSS_4h_fem_hyper, CpG_anno_4h_fem, "MF")
DMS_24h_fem_gostat_MF_hyper <- run.gostat(DMS_TSS_24h_fem_hyper, CpG_anno_24h_fem, "MF")
DMS_72h_fem_gostat_MF_hyper <- run.gostat(DMS_TSS_72h_fem_hyper, CpG_anno_72h_fem, "MF")

DMS_all_mal_gostat_MF_hyper <- run.gostat(DMS_TSS_all_mal_hyper, CpG_anno_all_mal, "MF")
DMS_05h_mal_gostat_MF_hyper <- run.gostat(DMS_TSS_05h_mal_hyper, CpG_anno_05h_mal, "MF")
DMS_1h_mal_gostat_MF_hyper <- run.gostat(DMS_TSS_1h_mal_hyper, CpG_anno_1h_mal, "MF")
DMS_4h_mal_gostat_MF_hyper <- run.gostat(DMS_TSS_4h_mal_hyper, CpG_anno_4h_mal, "MF")
DMS_24h_mal_gostat_MF_hyper <- run.gostat(DMS_TSS_24h_mal_hyper, CpG_anno_24h_mal, "MF")
DMS_72h_mal_gostat_MF_hyper <- run.gostat(DMS_TSS_72h_mal_hyper, CpG_anno_72h_mal, "MF")

DMS_all_fem_gostat_CC_hyper <- run.gostat(DMS_TSS_all_fem_hyper, CpG_anno_all_fem, "CC")
DMS_05h_fem_gostat_CC_hyper <- run.gostat(DMS_TSS_05h_fem_hyper, CpG_anno_05h_fem, "CC")
DMS_1h_fem_gostat_CC_hyper <- run.gostat(DMS_TSS_1h_fem_hyper, CpG_anno_1h_fem, "CC")
DMS_4h_fem_gostat_CC_hyper<- run.gostat(DMS_TSS_4h_fem_hyper, CpG_anno_4h_fem, "CC")
DMS_24h_fem_gostat_CC_hyper <- run.gostat(DMS_TSS_24h_fem_hyper, CpG_anno_24h_fem, "CC")
DMS_72h_fem_gostat_CC_hyper <- run.gostat(DMS_TSS_72h_fem_hyper, CpG_anno_72h_fem, "CC")

DMS_all_mal_gostat_CC_hyper <- run.gostat(DMS_TSS_all_mal_hyper, CpG_anno_all_mal, "CC")
DMS_05h_mal_gostat_CC_hyper <- run.gostat(DMS_TSS_05h_mal_hyper, CpG_anno_05h_mal, "CC")
DMS_1h_mal_gostat_CC_hyper <- run.gostat(DMS_TSS_1h_mal_hyper, CpG_anno_1h_mal, "CC")
DMS_4h_mal_gostat_CC_hyper <- run.gostat(DMS_TSS_4h_mal_hyper, CpG_anno_4h_mal, "CC")
DMS_24h_mal_gostat_CC_hyper <- run.gostat(DMS_TSS_24h_mal_hyper, CpG_anno_24h_mal, "CC")
DMS_72h_mal_gostat_CC_hyper <- run.gostat(DMS_TSS_72h_mal_hyper, CpG_anno_72h_mal, "CC")

DMR_all_fem_gostat_BP_hyper <- run.gostat(DMR_TSS_all_fem_hyper, regions_anno_all_fem, "BP")
DMR_05h_fem_gostat_BP_hyper <- run.gostat(DMR_TSS_05h_fem_hyper, regions_anno_05h_fem, "BP")
DMR_1h_fem_gostat_BP_hyper <- run.gostat(DMR_TSS_1h_fem_hyper, regions_anno_1h_fem, "BP")
DMR_4h_fem_gostat_BP_hyper<- run.gostat(DMR_TSS_4h_fem_hyper, regions_anno_4h_fem, "BP")
DMR_24h_fem_gostat_BP_hyper <- run.gostat(DMR_TSS_24h_fem_hyper, regions_anno_24h_fem, "BP")
DMR_72h_fem_gostat_BP_hyper <- run.gostat(DMR_TSS_72h_fem_hyper, regions_anno_72h_fem, "BP")

DMR_all_mal_gostat_BP_hyper <- run.gostat(DMR_TSS_all_mal_hyper, regions_anno_all_mal, "BP")
DMR_05h_mal_gostat_BP_hyper <- run.gostat(DMR_TSS_05h_mal_hyper, regions_anno_05h_mal, "BP")
DMR_1h_mal_gostat_BP_hyper <- run.gostat(DMR_TSS_1h_mal_hyper, regions_anno_1h_mal, "BP")
DMR_4h_mal_gostat_BP_hyper <- run.gostat(DMR_TSS_4h_mal_hyper, regions_anno_4h_mal, "BP")
DMR_24h_mal_gostat_BP_hyper <- run.gostat(DMR_TSS_24h_mal_hyper, regions_anno_24h_mal, "BP")
DMR_72h_mal_gostat_BP_hyper <- run.gostat(DMR_TSS_72h_mal_hyper, regions_anno_72h_mal, "BP")

DMR_all_fem_gostat_MF_hyper <- run.gostat(DMR_TSS_all_fem_hyper, regions_anno_all_fem, "MF")
DMR_05h_fem_gostat_MF_hyper <- run.gostat(DMR_TSS_05h_fem_hyper, regions_anno_05h_fem, "MF")
DMR_1h_fem_gostat_MF_hyper <- run.gostat(DMR_TSS_1h_fem_hyper, regions_anno_1h_fem, "MF")
DMR_4h_fem_gostat_MF_hyper<- run.gostat(DMR_TSS_4h_fem_hyper, regions_anno_4h_fem, "MF")
DMR_24h_fem_gostat_MF_hyper <- run.gostat(DMR_TSS_24h_fem_hyper, regions_anno_24h_fem, "MF")
DMR_72h_fem_gostat_MF_hyper <- run.gostat(DMR_TSS_72h_fem_hyper, regions_anno_72h_fem, "MF")

DMR_all_mal_gostat_MF_hyper <- run.gostat(DMR_TSS_all_mal_hyper, regions_anno_all_mal, "MF")
DMR_05h_mal_gostat_MF_hyper <- run.gostat(DMR_TSS_05h_mal_hyper, regions_anno_05h_mal, "MF")
DMR_1h_mal_gostat_MF_hyper <- run.gostat(DMR_TSS_1h_mal_hyper, regions_anno_1h_mal, "MF")
DMR_4h_mal_gostat_MF_hyper <- run.gostat(DMR_TSS_4h_mal_hyper, regions_anno_4h_mal, "MF")
DMR_24h_mal_gostat_MF_hyper <- run.gostat(DMR_TSS_24h_mal_hyper, regions_anno_24h_mal, "MF")
DMR_72h_mal_gostat_MF_hyper <- run.gostat(DMR_TSS_72h_mal_hyper, regions_anno_72h_mal, "MF")

DMR_all_fem_gostat_CC_hyper <- run.gostat(DMR_TSS_all_fem_hyper, regions_anno_all_fem, "CC")
DMR_05h_fem_gostat_CC_hyper <- run.gostat(DMR_TSS_05h_fem_hyper, regions_anno_05h_fem, "CC")
DMR_1h_fem_gostat_CC_hyper <- run.gostat(DMR_TSS_1h_fem_hyper, regions_anno_1h_fem, "CC")
DMR_4h_fem_gostat_CC_hyper<- run.gostat(DMR_TSS_4h_fem_hyper, regions_anno_4h_fem, "CC")
DMR_24h_fem_gostat_CC_hyper <- run.gostat(DMR_TSS_24h_fem_hyper, regions_anno_24h_fem, "CC")
DMR_72h_fem_gostat_CC_hyper <- run.gostat(DMR_TSS_72h_fem_hyper, regions_anno_72h_fem, "CC")

DMR_all_mal_gostat_CC_hyper <- run.gostat(DMR_TSS_all_mal_hyper, regions_anno_all_mal, "CC")
DMR_05h_mal_gostat_CC_hyper <- run.gostat(DMR_TSS_05h_mal_hyper, regions_anno_05h_mal, "CC")
DMR_1h_mal_gostat_CC_hyper <- run.gostat(DMR_TSS_1h_mal_hyper, regions_anno_1h_mal, "CC")
DMR_4h_mal_gostat_CC_hyper <- run.gostat(DMR_TSS_4h_mal_hyper, regions_anno_4h_mal, "CC")
DMR_24h_mal_gostat_CC_hyper <- run.gostat(DMR_TSS_24h_mal_hyper, regions_anno_24h_mal, "CC")
DMR_72h_mal_gostat_CC_hyper <- run.gostat(DMR_TSS_72h_mal_hyper, regions_anno_72h_mal, "CC")

#run on hypo gene sets 
DMS_all_fem_gostat_BP_hypo <- run.gostat(DMS_TSS_all_fem_hypo, CpG_anno_all_fem, "BP")
DMS_05h_fem_gostat_BP_hypo <- run.gostat(DMS_TSS_05h_fem_hypo, CpG_anno_05h_fem, "BP")
DMS_1h_fem_gostat_BP_hypo <- run.gostat(DMS_TSS_1h_fem_hypo, CpG_anno_1h_fem, "BP")
DMS_4h_fem_gostat_BP_hypo<- run.gostat(DMS_TSS_4h_fem_hypo, CpG_anno_4h_fem, "BP")
DMS_24h_fem_gostat_BP_hypo <- run.gostat(DMS_TSS_24h_fem_hypo, CpG_anno_24h_fem, "BP")
DMS_72h_fem_gostat_BP_hypo <- run.gostat(DMS_TSS_72h_fem_hypo, CpG_anno_72h_fem, "BP")

DMS_all_mal_gostat_BP_hypo <- run.gostat(DMS_TSS_all_mal_hypo, CpG_anno_all_mal, "BP")
DMS_05h_mal_gostat_BP_hypo <- run.gostat(DMS_TSS_05h_mal_hypo, CpG_anno_05h_mal, "BP")
DMS_1h_mal_gostat_BP_hypo <- run.gostat(DMS_TSS_1h_mal_hypo, CpG_anno_1h_mal, "BP")
DMS_4h_mal_gostat_BP_hypo <- run.gostat(DMS_TSS_4h_mal_hypo, CpG_anno_4h_mal, "BP")
DMS_24h_mal_gostat_BP_hypo <- run.gostat(DMS_TSS_24h_mal_hypo, CpG_anno_24h_mal, "BP")
DMS_72h_mal_gostat_BP_hypo <- run.gostat(DMS_TSS_72h_mal_hypo, CpG_anno_72h_mal, "BP")

DMS_all_fem_gostat_MF_hypo <- run.gostat(DMS_TSS_all_fem_hypo, CpG_anno_all_fem, "MF")
DMS_05h_fem_gostat_MF_hypo <- run.gostat(DMS_TSS_05h_fem_hypo, CpG_anno_05h_fem, "MF")
DMS_1h_fem_gostat_MF_hypo <- run.gostat(DMS_TSS_1h_fem_hypo, CpG_anno_1h_fem, "MF")
DMS_4h_fem_gostat_MF_hypo<- run.gostat(DMS_TSS_4h_fem_hypo, CpG_anno_4h_fem, "MF")
DMS_24h_fem_gostat_MF_hypo <- run.gostat(DMS_TSS_24h_fem_hypo, CpG_anno_24h_fem, "MF")
DMS_72h_fem_gostat_MF_hypo <- run.gostat(DMS_TSS_72h_fem_hypo, CpG_anno_72h_fem, "MF")

DMS_all_mal_gostat_MF_hypo <- run.gostat(DMS_TSS_all_mal_hypo, CpG_anno_all_mal, "MF")
DMS_05h_mal_gostat_MF_hypo <- run.gostat(DMS_TSS_05h_mal_hypo, CpG_anno_05h_mal, "MF")
DMS_1h_mal_gostat_MF_hypo <- run.gostat(DMS_TSS_1h_mal_hypo, CpG_anno_1h_mal, "MF")
DMS_4h_mal_gostat_MF_hypo <- run.gostat(DMS_TSS_4h_mal_hypo, CpG_anno_4h_mal, "MF")
DMS_24h_mal_gostat_MF_hypo <- run.gostat(DMS_TSS_24h_mal_hypo, CpG_anno_24h_mal, "MF")
DMS_72h_mal_gostat_MF_hypo <- run.gostat(DMS_TSS_72h_mal_hypo, CpG_anno_72h_mal, "MF")

DMS_all_fem_gostat_CC_hypo <- run.gostat(DMS_TSS_all_fem_hypo, CpG_anno_all_fem, "CC")
DMS_05h_fem_gostat_CC_hypo <- run.gostat(DMS_TSS_05h_fem_hypo, CpG_anno_05h_fem, "CC")
DMS_1h_fem_gostat_CC_hypo <- run.gostat(DMS_TSS_1h_fem_hypo, CpG_anno_1h_fem, "CC")
DMS_4h_fem_gostat_CC_hypo<- run.gostat(DMS_TSS_4h_fem_hypo, CpG_anno_4h_fem, "CC")
DMS_24h_fem_gostat_CC_hypo <- run.gostat(DMS_TSS_24h_fem_hypo, CpG_anno_24h_fem, "CC")
DMS_72h_fem_gostat_CC_hypo <- run.gostat(DMS_TSS_72h_fem_hypo, CpG_anno_72h_fem, "CC")

DMS_all_mal_gostat_CC_hypo <- run.gostat(DMS_TSS_all_mal_hypo, CpG_anno_all_mal, "CC")
DMS_05h_mal_gostat_CC_hypo <- run.gostat(DMS_TSS_05h_mal_hypo, CpG_anno_05h_mal, "CC")
DMS_1h_mal_gostat_CC_hypo <- run.gostat(DMS_TSS_1h_mal_hypo, CpG_anno_1h_mal, "CC")
DMS_4h_mal_gostat_CC_hypo <- run.gostat(DMS_TSS_4h_mal_hypo, CpG_anno_4h_mal, "CC")
DMS_24h_mal_gostat_CC_hypo <- run.gostat(DMS_TSS_24h_mal_hypo, CpG_anno_24h_mal, "CC")
DMS_72h_mal_gostat_CC_hypo <- run.gostat(DMS_TSS_72h_mal_hypo, CpG_anno_72h_mal, "CC")

DMR_all_fem_gostat_BP_hypo <- run.gostat(DMR_TSS_all_fem_hypo, regions_anno_all_fem, "BP")
DMR_05h_fem_gostat_BP_hypo <- run.gostat(DMR_TSS_05h_fem_hypo, regions_anno_05h_fem, "BP")
DMR_1h_fem_gostat_BP_hypo <- run.gostat(DMR_TSS_1h_fem_hypo, regions_anno_1h_fem, "BP")
DMR_4h_fem_gostat_BP_hypo<- run.gostat(DMR_TSS_4h_fem_hypo, regions_anno_4h_fem, "BP")
DMR_24h_fem_gostat_BP_hypo <- run.gostat(DMR_TSS_24h_fem_hypo, regions_anno_24h_fem, "BP")
DMR_72h_fem_gostat_BP_hypo <- run.gostat(DMR_TSS_72h_fem_hypo, regions_anno_72h_fem, "BP")

DMR_all_mal_gostat_BP_hypo <- run.gostat(DMR_TSS_all_mal_hypo, regions_anno_all_mal, "BP")
DMR_05h_mal_gostat_BP_hypo <- run.gostat(DMR_TSS_05h_mal_hypo, regions_anno_05h_mal, "BP")
DMR_1h_mal_gostat_BP_hypo <- run.gostat(DMR_TSS_1h_mal_hypo, regions_anno_1h_mal, "BP")
DMR_4h_mal_gostat_BP_hypo <- run.gostat(DMR_TSS_4h_mal_hypo, regions_anno_4h_mal, "BP")
DMR_24h_mal_gostat_BP_hypo <- run.gostat(DMR_TSS_24h_mal_hypo, regions_anno_24h_mal, "BP")
DMR_72h_mal_gostat_BP_hypo <- run.gostat(DMR_TSS_72h_mal_hypo, regions_anno_72h_mal, "BP")

DMR_all_fem_gostat_MF_hypo <- run.gostat(DMR_TSS_all_fem_hypo, regions_anno_all_fem, "MF")
DMR_05h_fem_gostat_MF_hypo <- run.gostat(DMR_TSS_05h_fem_hypo, regions_anno_05h_fem, "MF")
DMR_1h_fem_gostat_MF_hypo <- run.gostat(DMR_TSS_1h_fem_hypo, regions_anno_1h_fem, "MF")
DMR_4h_fem_gostat_MF_hypo<- run.gostat(DMR_TSS_4h_fem_hypo, regions_anno_4h_fem, "MF")
DMR_24h_fem_gostat_MF_hypo <- run.gostat(DMR_TSS_24h_fem_hypo, regions_anno_24h_fem, "MF")
DMR_72h_fem_gostat_MF_hypo <- run.gostat(DMR_TSS_72h_fem_hypo, regions_anno_72h_fem, "MF")

DMR_all_mal_gostat_MF_hypo <- run.gostat(DMR_TSS_all_mal_hypo, regions_anno_all_mal, "MF")
DMR_05h_mal_gostat_MF_hypo <- run.gostat(DMR_TSS_05h_mal_hypo, regions_anno_05h_mal, "MF")
DMR_1h_mal_gostat_MF_hypo <- run.gostat(DMR_TSS_1h_mal_hypo, regions_anno_1h_mal, "MF")
DMR_4h_mal_gostat_MF_hypo <- run.gostat(DMR_TSS_4h_mal_hypo, regions_anno_4h_mal, "MF")
DMR_24h_mal_gostat_MF_hypo <- run.gostat(DMR_TSS_24h_mal_hypo, regions_anno_24h_mal, "MF")
DMR_72h_mal_gostat_MF_hypo <- run.gostat(DMR_TSS_72h_mal_hypo, regions_anno_72h_mal, "MF")

DMR_all_fem_gostat_CC_hypo <- run.gostat(DMR_TSS_all_fem_hypo, regions_anno_all_fem, "CC")
DMR_05h_fem_gostat_CC_hypo <- run.gostat(DMR_TSS_05h_fem_hypo, regions_anno_05h_fem, "CC")
DMR_1h_fem_gostat_CC_hypo <- run.gostat(DMR_TSS_1h_fem_hypo, regions_anno_1h_fem, "CC")
DMR_4h_fem_gostat_CC_hypo<- run.gostat(DMR_TSS_4h_fem_hypo, regions_anno_4h_fem, "CC")
DMR_24h_fem_gostat_CC_hypo <- run.gostat(DMR_TSS_24h_fem_hypo, regions_anno_24h_fem, "CC")
DMR_72h_fem_gostat_CC_hypo <- run.gostat(DMR_TSS_72h_fem_hypo, regions_anno_72h_fem, "CC")

DMR_all_mal_gostat_CC_hypo <- run.gostat(DMR_TSS_all_mal_hypo, regions_anno_all_mal, "CC")
DMR_05h_mal_gostat_CC_hypo <- run.gostat(DMR_TSS_05h_mal_hypo, regions_anno_05h_mal, "CC")
DMR_1h_mal_gostat_CC_hypo <- run.gostat(DMR_TSS_1h_mal_hypo, regions_anno_1h_mal, "CC")
DMR_4h_mal_gostat_CC_hypo <- run.gostat(DMR_TSS_4h_mal_hypo, regions_anno_4h_mal, "CC")
DMR_24h_mal_gostat_CC_hypo <- run.gostat(DMR_TSS_24h_mal_hypo, regions_anno_24h_mal, "CC")
DMR_72h_mal_gostat_CC_hypo <- run.gostat(DMR_TSS_72h_mal_hypo, regions_anno_72h_mal, "CC")

#save results as a table
write.csv(DMS_all_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_all_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_05h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_1h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_4h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_24h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_72h_fem_gostat_BP.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_all_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_05h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_1h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_4h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_24h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMS_72h_mal_gostat_BP.csv", row.names = FALSE)

write.csv(DMS_all_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_all_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_05h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_1h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_4h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_24h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_72h_fem_gostat_MF.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_all_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_05h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_1h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_4h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_24h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMS_72h_mal_gostat_MF.csv", row.names = FALSE)

write.csv(DMS_all_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_all_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_05h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_1h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_4h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_24h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_72h_fem_gostat_CC.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_all_mal_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_05h_mal_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_1h_mal_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_4h_mal_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_24h_mal_gostat_CC.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_CC, "./shortTerm_exp/data/go_res_od/DMS_72h_mal_gostat_CC.csv", row.names = FALSE)

write.csv(DMR_all_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_all_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_05h_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_05h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_1h_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_1h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_4h_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_4h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_24h_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_24h_fem_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_72h_fem_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_72h_fem_gostat_BP.csv", row.names = FALSE)

write.csv(DMR_all_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_all_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_05h_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_05h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_1h_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_1h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_4h_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_4h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_24h_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_24h_mal_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_72h_mal_gostat_BP, "./shortTerm_exp/data/go_res_od/DMR_72h_mal_gostat_BP.csv", row.names = FALSE)

write.csv(DMR_all_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_all_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_05h_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_05h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_1h_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_1h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_4h_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_4h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_24h_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_24h_fem_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_72h_fem_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_72h_fem_gostat_MF.csv", row.names = FALSE)

write.csv(DMR_all_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_all_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_05h_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_05h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_1h_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_1h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_4h_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_4h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_24h_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_24h_mal_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_72h_mal_gostat_MF, "./shortTerm_exp/data/go_res_od/DMR_72h_mal_gostat_MF.csv", row.names = FALSE)

write.csv(DMR_all_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMR_all_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMR_05h_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMR_05h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMR_1h_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMR_1h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMR_4h_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMR_4h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMR_24h_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMR_24h_fem_gostat_CC.csv", row.names = FALSE)
write.csv(DMR_72h_fem_gostat_CC, "./shortTerm_exp/data/go_res_od/DMR_72h_fem_gostat_CC.csv", row.names = FALSE)

#hyper
write.csv(DMS_all_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_all_fem_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_05h_fem_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_1h_fem_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_4h_fem_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_24h_fem_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_72h_fem_gostat_BP_hyper.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_all_mal_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_05h_mal_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_1h_mal_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_4h_mal_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_24h_mal_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_72h_mal_gostat_BP_hyper.csv", row.names = FALSE)

write.csv(DMS_all_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_all_fem_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_05h_fem_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_1h_fem_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_4h_fem_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_24h_fem_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_72h_fem_gostat_MF_hyper.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_all_mal_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_05h_mal_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_1h_mal_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_4h_mal_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_24h_mal_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_72h_mal_gostat_MF_hyper.csv", row.names = FALSE)

write.csv(DMS_all_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_all_fem_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_05h_fem_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_1h_fem_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_4h_fem_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_24h_fem_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_72h_fem_gostat_CC_hyper.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_all_mal_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_05h_mal_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_1h_mal_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_4h_mal_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_24h_mal_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_72h_mal_gostat_CC_hyper.csv", row.names = FALSE)

write.csv(DMR_all_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_all_fem_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_05h_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_05h_fem_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_1h_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_1h_fem_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_4h_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_4h_fem_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_24h_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_24h_fem_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_72h_fem_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_72h_fem_gostat_BP_hyper.csv", row.names = FALSE)

write.csv(DMR_all_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_all_mal_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_05h_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_05h_mal_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_1h_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_1h_mal_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_4h_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_4h_mal_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_24h_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_24h_mal_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_72h_mal_gostat_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_72h_mal_gostat_BP_hyper.csv", row.names = FALSE)

write.csv(DMR_all_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_all_fem_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_05h_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_05h_fem_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_1h_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_1h_fem_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_4h_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_4h_fem_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_24h_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_24h_fem_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_72h_fem_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_72h_fem_gostat_MF_hyper.csv", row.names = FALSE)

write.csv(DMR_all_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_all_mal_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_05h_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_05h_mal_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_1h_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_1h_mal_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_4h_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_4h_mal_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_24h_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_24h_mal_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_72h_mal_gostat_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_72h_mal_gostat_MF_hyper.csv", row.names = FALSE)

write.csv(DMR_all_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_all_fem_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMR_05h_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_05h_fem_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMR_1h_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_1h_fem_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMR_4h_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_4h_fem_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMR_24h_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_24h_fem_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMR_72h_fem_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_72h_fem_gostat_CC_hyper.csv", row.names = FALSE)

write.csv(DMR_all_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_all_mal_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMR_05h_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_05h_mal_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMR_1h_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_1h_mal_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMR_4h_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_4h_mal_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMR_24h_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_24h_mal_gostat_CC_hyper.csv", row.names = FALSE)
write.csv(DMR_72h_mal_gostat_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_72h_mal_gostat_CC_hyper.csv", row.names = FALSE)

#hypo
write.csv(DMS_all_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_all_fem_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_05h_fem_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_1h_fem_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_4h_fem_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_24h_fem_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_72h_fem_gostat_BP_hypo.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_all_mal_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_05h_mal_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_1h_mal_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_4h_mal_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_24h_mal_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_72h_mal_gostat_BP_hypo.csv", row.names = FALSE)

write.csv(DMS_all_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_all_fem_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_05h_fem_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_1h_fem_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_4h_fem_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_24h_fem_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_72h_fem_gostat_MF_hypo.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_all_mal_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_05h_mal_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_1h_mal_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_4h_mal_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_24h_mal_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_72h_mal_gostat_MF_hypo.csv", row.names = FALSE)

write.csv(DMS_all_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_all_fem_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMS_05h_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_05h_fem_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMS_1h_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_1h_fem_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMS_4h_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_4h_fem_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMS_24h_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_24h_fem_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMS_72h_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_72h_fem_gostat_CC_hypo.csv", row.names = FALSE)

write.csv(DMS_all_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_all_mal_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMS_05h_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_05h_mal_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMS_1h_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_1h_mal_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMS_4h_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_4h_mal_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMS_24h_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_24h_mal_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMS_72h_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_72h_mal_gostat_CC_hypo.csv", row.names = FALSE)

write.csv(DMR_all_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_all_fem_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_05h_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_05h_fem_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_1h_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_1h_fem_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_4h_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_4h_fem_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_24h_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_24h_fem_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_72h_fem_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_72h_fem_gostat_BP_hypo.csv", row.names = FALSE)

write.csv(DMR_all_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_all_mal_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_05h_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_05h_mal_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_1h_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_1h_mal_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_4h_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_4h_mal_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_24h_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_24h_mal_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_72h_mal_gostat_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_72h_mal_gostat_BP_hypo.csv", row.names = FALSE)

write.csv(DMR_all_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_all_fem_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_05h_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_05h_fem_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_1h_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_1h_fem_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_4h_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_4h_fem_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_24h_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_24h_fem_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_72h_fem_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_72h_fem_gostat_MF_hypo.csv", row.names = FALSE)

write.csv(DMR_all_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_all_mal_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_05h_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_05h_mal_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_1h_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_1h_mal_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_4h_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_4h_mal_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_24h_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_24h_mal_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_72h_mal_gostat_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_72h_mal_gostat_MF_hypo.csv", row.names = FALSE)

write.csv(DMR_all_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_all_fem_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMR_05h_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_05h_fem_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMR_1h_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_1h_fem_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMR_4h_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_4h_fem_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMR_24h_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_24h_fem_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMR_72h_fem_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_72h_fem_gostat_CC_hypo.csv", row.names = FALSE)

write.csv(DMR_all_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_all_mal_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMR_05h_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_05h_mal_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMR_1h_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_1h_mal_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMR_4h_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_4h_mal_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMR_24h_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_24h_mal_gostat_CC_hypo.csv", row.names = FALSE)
write.csv(DMR_72h_mal_gostat_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_72h_mal_gostat_CC_hypo.csv", row.names = FALSE)

#merge dataframse and save
DMS_05h_fem_gostat_BP$time_point <- "0.5h"
DMS_1h_fem_gostat_BP$time_point <- "1h"
DMS_4h_fem_gostat_BP$time_point <- "4h"
DMS_24h_fem_gostat_BP$time_point <- "24h"
DMS_72h_fem_gostat_BP$time_point <- "72h"

DMS_05h_fem_gostat_MF$time_point <- "0.5h"
DMS_1h_fem_gostat_MF$time_point <- "1h"
DMS_4h_fem_gostat_MF$time_point <- "4h"
DMS_24h_fem_gostat_MF$time_point <- "24h"
DMS_72h_fem_gostat_MF$time_point <- "72h"

DMS_05h_fem_gostat_CC$time_point <- "0.5h"
DMS_1h_fem_gostat_CC$time_point <- "1h"
DMS_4h_fem_gostat_CC$time_point <- "4h"
DMS_24h_fem_gostat_CC$time_point <- "24h"
DMS_72h_fem_gostat_CC$time_point <- "72h"

DMS_fem_all_BP <- rbind(DMS_05h_fem_gostat_BP,DMS_1h_fem_gostat_BP,DMS_4h_fem_gostat_BP,DMS_24h_fem_gostat_BP,DMS_72h_fem_gostat_BP)
DMS_fem_all_MF <- rbind(DMS_05h_fem_gostat_MF,DMS_1h_fem_gostat_MF,DMS_4h_fem_gostat_MF,DMS_24h_fem_gostat_MF,DMS_72h_fem_gostat_MF)
DMS_fem_all_CC <- rbind(DMS_05h_fem_gostat_CC,DMS_1h_fem_gostat_CC,DMS_4h_fem_gostat_CC,DMS_24h_fem_gostat_CC,DMS_72h_fem_gostat_CC)

write.csv(DMS_fem_all_BP, "./shortTerm_exp/data/go_res_od/DMS_fem_combined_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_fem_all_MF, "./shortTerm_exp/data/go_res_od/DMS_fem_combined_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_fem_all_CC, "./shortTerm_exp/data/go_res_od/DMS_fem_combined_gostat_CC.csv", row.names = FALSE)

DMS_05h_mal_gostat_BP$time_point <- "0.5h"
DMS_1h_mal_gostat_BP$time_point <- "1h"
DMS_4h_mal_gostat_BP$time_point <- "4h"
DMS_24h_mal_gostat_BP$time_point <- "24h"
DMS_72h_mal_gostat_BP$time_point <- "72h"

DMS_05h_mal_gostat_MF$time_point <- "0.5h"
DMS_1h_mal_gostat_MF$time_point <- "1h"
DMS_4h_mal_gostat_MF$time_point <- "4h"
DMS_24h_mal_gostat_MF$time_point <- "24h"
DMS_72h_mal_gostat_MF$time_point <- "72h"

DMS_05h_mal_gostat_CC$time_point <- "0.5h"
DMS_1h_mal_gostat_CC$time_point <- "1h"
DMS_4h_mal_gostat_CC$time_point <- "4h"
DMS_24h_mal_gostat_CC$time_point <- "24h"
DMS_72h_mal_gostat_CC$time_point <- "72h"

DMS_mal_all_BP <- rbind(DMS_05h_mal_gostat_BP,DMS_1h_mal_gostat_BP,DMS_4h_mal_gostat_BP,DMS_24h_mal_gostat_BP,DMS_72h_mal_gostat_BP)
DMS_mal_all_MF <- rbind(DMS_05h_mal_gostat_MF,DMS_1h_mal_gostat_MF,DMS_4h_mal_gostat_MF,DMS_24h_mal_gostat_MF,DMS_72h_mal_gostat_MF)
DMS_mal_all_CC <- rbind(DMS_05h_mal_gostat_CC,DMS_1h_mal_gostat_CC,DMS_4h_mal_gostat_CC,DMS_24h_mal_gostat_CC,DMS_72h_mal_gostat_CC)

write.csv(DMS_mal_all_BP, "./shortTerm_exp/data/go_res_od/DMS_mal_combined_gostat_BP.csv", row.names = FALSE)
write.csv(DMS_mal_all_MF, "./shortTerm_exp/data/go_res_od/DMS_mal_combined_gostat_MF.csv", row.names = FALSE)
write.csv(DMS_mal_all_CC, "./shortTerm_exp/data/go_res_od/DMS_mal_combined_gostat_CC.csv", row.names = FALSE)

DMR_05h_fem_gostat_BP$time_point <- "0.5h"
DMR_1h_fem_gostat_BP$time_point <- "1h"
DMR_4h_fem_gostat_BP$time_point <- "4h"
DMR_24h_fem_gostat_BP$time_point <- "24h"
DMR_72h_fem_gostat_BP$time_point <- "72h"

DMR_05h_fem_gostat_MF$time_point <- "0.5h"
DMR_1h_fem_gostat_MF$time_point <- "1h"
DMR_4h_fem_gostat_MF$time_point <- "4h"
DMR_24h_fem_gostat_MF$time_point <- "24h"
DMR_72h_fem_gostat_MF$time_point <- "72h"

DMR_05h_fem_gostat_CC$time_point <- "0.5h"
DMR_1h_fem_gostat_CC$time_point <- "1h"
DMR_4h_fem_gostat_CC$time_point <- "4h"
DMR_24h_fem_gostat_CC$time_point <- "24h"
DMR_72h_fem_gostat_CC$time_point <- "72h"

DMR_fem_all_BP <- rbind(DMR_05h_fem_gostat_BP,DMR_1h_fem_gostat_BP,DMR_4h_fem_gostat_BP,DMR_24h_fem_gostat_BP,DMR_72h_fem_gostat_BP)
DMR_fem_all_MF <- rbind(DMR_05h_fem_gostat_MF,DMR_1h_fem_gostat_MF,DMR_4h_fem_gostat_MF,DMR_24h_fem_gostat_MF,DMR_72h_fem_gostat_MF)
DMR_fem_all_CC <- rbind(DMR_05h_fem_gostat_CC,DMR_1h_fem_gostat_CC,DMR_4h_fem_gostat_CC,DMR_24h_fem_gostat_CC,DMR_72h_fem_gostat_CC)

write.csv(DMR_fem_all_BP, "./shortTerm_exp/data/go_res_od/DMR_fem_combined_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_fem_all_MF, "./shortTerm_exp/data/go_res_od/DMR_fem_combined_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_fem_all_CC, "./shortTerm_exp/data/go_res_od/DMR_fem_combined_gostat_CC.csv", row.names = FALSE)

DMR_05h_mal_gostat_BP$time_point <- "0.5h"
DMR_1h_mal_gostat_BP$time_point <- "1h"
DMR_4h_mal_gostat_BP$time_point <- "4h"
DMR_24h_mal_gostat_BP$time_point <- "24h"
DMR_72h_mal_gostat_BP$time_point <- "72h"

DMR_05h_mal_gostat_MF$time_point <- "0.5h"
DMR_1h_mal_gostat_MF$time_point <- "1h"
DMR_4h_mal_gostat_MF$time_point <- "4h"
DMR_24h_mal_gostat_MF$time_point <- "24h"
DMR_72h_mal_gostat_MF$time_point <- "72h"

DMR_05h_mal_gostat_CC$time_point <- "0.5h"
DMR_1h_mal_gostat_CC$time_point <- "1h"
DMR_4h_mal_gostat_CC$time_point <- "4h"
DMR_24h_mal_gostat_CC$time_point <- "24h"
DMR_72h_mal_gostat_CC$time_point <- "72h"

DMR_mal_all_BP <- rbind(DMR_05h_mal_gostat_BP,DMR_1h_mal_gostat_BP,DMR_4h_mal_gostat_BP,DMR_24h_mal_gostat_BP,DMR_72h_mal_gostat_BP)
DMR_mal_all_MF <- rbind(DMR_05h_mal_gostat_MF,DMR_1h_mal_gostat_MF,DMR_4h_mal_gostat_MF,DMR_24h_mal_gostat_MF,DMR_72h_mal_gostat_MF)
DMR_mal_all_CC <- rbind(DMR_05h_mal_gostat_CC,DMR_1h_mal_gostat_CC,DMR_4h_mal_gostat_CC,DMR_24h_mal_gostat_CC,DMR_72h_mal_gostat_CC)

write.csv(DMR_mal_all_BP, "./shortTerm_exp/data/go_res_od/DMR_mal_combined_gostat_BP.csv", row.names = FALSE)
write.csv(DMR_mal_all_MF, "./shortTerm_exp/data/go_res_od/DMR_mal_combined_gostat_MF.csv", row.names = FALSE)
write.csv(DMR_mal_all_CC, "./shortTerm_exp/data/go_res_od/DMR_mal_combined_gostat_CC.csv", row.names = FALSE)

#hyper
DMS_05h_fem_gostat_BP_hyper$time_point <- "0.5h"
DMS_1h_fem_gostat_BP_hyper$time_point <- "1h"
DMS_4h_fem_gostat_BP_hyper$time_point <- "4h"
DMS_24h_fem_gostat_BP_hyper$time_point <- "24h"
DMS_72h_fem_gostat_BP_hyper$time_point <- "72h"

DMS_05h_fem_gostat_MF_hyper$time_point <- "0.5h"
DMS_1h_fem_gostat_MF_hyper$time_point <- "1h"
DMS_4h_fem_gostat_MF_hyper$time_point <- "4h"
DMS_24h_fem_gostat_MF_hyper$time_point <- "24h"
DMS_72h_fem_gostat_MF_hyper$time_point <- "72h"

DMS_05h_fem_gostat_CC_hyper$time_point <- "0.5h"
DMS_1h_fem_gostat_CC_hyper$time_point <- "1h"
DMS_4h_fem_gostat_CC_hyper$time_point <- "4h"
DMS_24h_fem_gostat_CC_hyper$time_point <- "24h"
DMS_72h_fem_gostat_CC_hyper$time_point <- "72h"

DMS_fem_all_BP_hyper <- rbind(DMS_05h_fem_gostat_BP_hyper,DMS_1h_fem_gostat_BP_hyper,DMS_4h_fem_gostat_BP_hyper,DMS_24h_fem_gostat_BP_hyper,DMS_72h_fem_gostat_BP_hyper)
DMS_fem_all_MF_hyper <- rbind(DMS_05h_fem_gostat_MF_hyper,DMS_1h_fem_gostat_MF_hyper,DMS_4h_fem_gostat_MF_hyper,DMS_24h_fem_gostat_MF_hyper,DMS_72h_fem_gostat_MF_hyper)
DMS_fem_all_CC_hyper <- rbind(DMS_05h_fem_gostat_CC_hyper,DMS_1h_fem_gostat_CC_hyper,DMS_4h_fem_gostat_CC_hyper,DMS_24h_fem_gostat_CC_hyper,DMS_72h_fem_gostat_CC_hyper)

write.csv(DMS_fem_all_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_fem_combined_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_fem_all_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_fem_combined_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_fem_all_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_fem_combined_gostat_CC_hyper.csv", row.names = FALSE)

DMS_05h_mal_gostat_BP_hyper$time_point <- "0.5h"
DMS_1h_mal_gostat_BP_hyper$time_point <- "1h"
DMS_4h_mal_gostat_BP_hyper$time_point <- "4h"
DMS_24h_mal_gostat_BP_hyper$time_point <- "24h"
DMS_72h_mal_gostat_BP_hyper$time_point <- "72h"

DMS_05h_mal_gostat_MF_hyper$time_point <- "0.5h"
DMS_1h_mal_gostat_MF_hyper$time_point <- "1h"
DMS_4h_mal_gostat_MF_hyper$time_point <- "4h"
DMS_24h_mal_gostat_MF_hyper$time_point <- "24h"
DMS_72h_mal_gostat_MF_hyper$time_point <- "72h"

DMS_05h_mal_gostat_CC_hyper$time_point <- "0.5h"
DMS_1h_mal_gostat_CC_hyper$time_point <- "1h"
DMS_4h_mal_gostat_CC_hyper$time_point <- "4h"
DMS_24h_mal_gostat_CC_hyper$time_point <- "24h"
DMS_72h_mal_gostat_CC_hyper$time_point <- "72h"

DMS_mal_all_BP_hyper <- rbind(DMS_05h_mal_gostat_BP_hyper,DMS_1h_mal_gostat_BP_hyper,DMS_4h_mal_gostat_BP_hyper,DMS_24h_mal_gostat_BP_hyper,DMS_72h_mal_gostat_BP_hyper)
DMS_mal_all_MF_hyper <- rbind(DMS_05h_mal_gostat_MF_hyper,DMS_1h_mal_gostat_MF_hyper,DMS_4h_mal_gostat_MF_hyper,DMS_24h_mal_gostat_MF_hyper,DMS_72h_mal_gostat_MF_hyper)
DMS_mal_all_CC_hyper <- rbind(DMS_05h_mal_gostat_CC_hyper,DMS_1h_mal_gostat_CC_hyper,DMS_4h_mal_gostat_CC_hyper,DMS_24h_mal_gostat_CC_hyper,DMS_72h_mal_gostat_CC_hyper)

write.csv(DMS_mal_all_BP_hyper, "./shortTerm_exp/data/go_res_od/DMS_mal_combined_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMS_mal_all_MF_hyper, "./shortTerm_exp/data/go_res_od/DMS_mal_combined_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMS_mal_all_CC_hyper, "./shortTerm_exp/data/go_res_od/DMS_mal_combined_gostat_CC_hyper.csv", row.names = FALSE)

DMR_05h_fem_gostat_BP_hyper$time_point <- "0.5h"
DMR_1h_fem_gostat_BP_hyper$time_point <- "1h"
DMR_4h_fem_gostat_BP_hyper$time_point <- "4h"
DMR_24h_fem_gostat_BP_hyper$time_point <- "24h"
DMR_72h_fem_gostat_BP_hyper$time_point <- "72h"

DMR_05h_fem_gostat_MF_hyper$time_point <- "0.5h"
DMR_1h_fem_gostat_MF_hyper$time_point <- "1h"
DMR_4h_fem_gostat_MF_hyper$time_point <- "4h"
DMR_24h_fem_gostat_MF_hyper$time_point <- "24h"
DMR_72h_fem_gostat_MF_hyper$time_point <- "72h"

DMR_05h_fem_gostat_CC_hyper$time_point <- "0.5h"
DMR_1h_fem_gostat_CC_hyper$time_point <- "1h"
DMR_4h_fem_gostat_CC_hyper$time_point <- "4h"
DMR_24h_fem_gostat_CC_hyper$time_point <- "24h"
DMR_72h_fem_gostat_CC_hyper$time_point <- "72h"

DMR_fem_all_BP_hyper <- rbind(DMR_05h_fem_gostat_BP_hyper,DMR_1h_fem_gostat_BP_hyper,DMR_4h_fem_gostat_BP_hyper,DMR_24h_fem_gostat_BP_hyper,DMR_72h_fem_gostat_BP_hyper)
DMR_fem_all_MF_hyper <- rbind(DMR_05h_fem_gostat_MF_hyper,DMR_1h_fem_gostat_MF_hyper,DMR_4h_fem_gostat_MF_hyper,DMR_24h_fem_gostat_MF_hyper,DMR_72h_fem_gostat_MF_hyper)
DMR_fem_all_CC_hyper <- rbind(DMR_05h_fem_gostat_CC_hyper,DMR_1h_fem_gostat_CC_hyper,DMR_4h_fem_gostat_CC_hyper,DMR_24h_fem_gostat_CC_hyper,DMR_72h_fem_gostat_CC_hyper)

write.csv(DMR_fem_all_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_fem_combined_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_fem_all_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_fem_combined_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_fem_all_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_fem_combined_gostat_CC_hyper.csv", row.names = FALSE)

DMR_05h_mal_gostat_BP_hyper$time_point <- "0.5h"
DMR_1h_mal_gostat_BP_hyper$time_point <- "1h"
DMR_4h_mal_gostat_BP_hyper$time_point <- "4h"
DMR_24h_mal_gostat_BP_hyper$time_point <- "24h"
DMR_72h_mal_gostat_BP_hyper$time_point <- "72h"

DMR_05h_mal_gostat_MF_hyper$time_point <- "0.5h"
DMR_1h_mal_gostat_MF_hyper$time_point <- "1h"
DMR_4h_mal_gostat_MF_hyper$time_point <- "4h"
DMR_24h_mal_gostat_MF_hyper$time_point <- "24h"
DMR_72h_mal_gostat_MF_hyper$time_point <- "72h"

DMR_05h_mal_gostat_CC_hyper$time_point <- "0.5h"
DMR_1h_mal_gostat_CC_hyper$time_point <- "1h"
DMR_4h_mal_gostat_CC_hyper$time_point <- "4h"
DMR_24h_mal_gostat_CC_hyper$time_point <- "24h"
DMR_72h_mal_gostat_CC_hyper$time_point <- "72h"

DMR_mal_all_BP_hyper <- rbind(DMR_05h_mal_gostat_BP_hyper,DMR_1h_mal_gostat_BP_hyper,DMR_4h_mal_gostat_BP_hyper,DMR_24h_mal_gostat_BP_hyper,DMR_72h_mal_gostat_BP_hyper)
DMR_mal_all_MF_hyper <- rbind(DMR_05h_mal_gostat_MF_hyper,DMR_1h_mal_gostat_MF_hyper,DMR_4h_mal_gostat_MF_hyper,DMR_24h_mal_gostat_MF_hyper,DMR_72h_mal_gostat_MF_hyper)
DMR_mal_all_CC_hyper <- rbind(DMR_05h_mal_gostat_CC_hyper,DMR_1h_mal_gostat_CC_hyper,DMR_4h_mal_gostat_CC_hyper,DMR_24h_mal_gostat_CC_hyper,DMR_72h_mal_gostat_CC_hyper)

write.csv(DMR_mal_all_BP_hyper, "./shortTerm_exp/data/go_res_od/DMR_mal_combined_gostat_BP_hyper.csv", row.names = FALSE)
write.csv(DMR_mal_all_MF_hyper, "./shortTerm_exp/data/go_res_od/DMR_mal_combined_gostat_MF_hyper.csv", row.names = FALSE)
write.csv(DMR_mal_all_CC_hyper, "./shortTerm_exp/data/go_res_od/DMR_mal_combined_gostat_CC_hyper.csv", row.names = FALSE)

#hypo
DMS_05h_fem_gostat_BP_hypo$time_point <- "0.5h"
DMS_1h_fem_gostat_BP_hypo$time_point <- "1h"
DMS_4h_fem_gostat_BP_hypo$time_point <- "4h"
DMS_24h_fem_gostat_BP_hypo$time_point <- "24h"
DMS_72h_fem_gostat_BP_hypo$time_point <- "72h"

DMS_05h_fem_gostat_MF_hypo$time_point <- "0.5h"
DMS_1h_fem_gostat_MF_hypo$time_point <- "1h"
DMS_4h_fem_gostat_MF_hypo$time_point <- "4h"
DMS_24h_fem_gostat_MF_hypo$time_point <- "24h"
DMS_72h_fem_gostat_MF_hypo$time_point <- "72h"

DMS_05h_fem_gostat_CC_hypo$time_point <- "0.5h"
DMS_1h_fem_gostat_CC_hypo$time_point <- "1h"
DMS_4h_fem_gostat_CC_hypo$time_point <- "4h"
DMS_24h_fem_gostat_CC_hypo$time_point <- "24h"
DMS_72h_fem_gostat_CC_hypo$time_point <- "72h"

DMS_fem_all_BP_hypo <- rbind(DMS_05h_fem_gostat_BP_hypo,DMS_1h_fem_gostat_BP_hypo,DMS_4h_fem_gostat_BP_hypo,DMS_24h_fem_gostat_BP_hypo,DMS_72h_fem_gostat_BP_hypo)
DMS_fem_all_MF_hypo <- rbind(DMS_05h_fem_gostat_MF_hypo,DMS_1h_fem_gostat_MF_hypo,DMS_4h_fem_gostat_MF_hypo,DMS_24h_fem_gostat_MF_hypo,DMS_72h_fem_gostat_MF_hypo)
DMS_fem_all_CC_hypo <- rbind(DMS_05h_fem_gostat_CC_hypo,DMS_1h_fem_gostat_CC_hypo,DMS_4h_fem_gostat_CC_hypo,DMS_24h_fem_gostat_CC_hypo,DMS_72h_fem_gostat_CC_hypo)

write.csv(DMS_fem_all_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_fem_combined_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_fem_all_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_fem_combined_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_fem_all_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_fem_combined_gostat_CC_hypo.csv", row.names = FALSE)

DMS_05h_mal_gostat_BP_hypo$time_point <- "0.5h"
DMS_1h_mal_gostat_BP_hypo$time_point <- "1h"
DMS_4h_mal_gostat_BP_hypo$time_point <- "4h"
DMS_24h_mal_gostat_BP_hypo$time_point <- "24h"
DMS_72h_mal_gostat_BP_hypo$time_point <- "72h"

DMS_05h_mal_gostat_MF_hypo$time_point <- "0.5h"
DMS_1h_mal_gostat_MF_hypo$time_point <- "1h"
DMS_4h_mal_gostat_MF_hypo$time_point <- "4h"
DMS_24h_mal_gostat_MF_hypo$time_point <- "24h"
DMS_72h_mal_gostat_MF_hypo$time_point <- "72h"

DMS_05h_mal_gostat_CC_hypo$time_point <- "0.5h"
DMS_1h_mal_gostat_CC_hypo$time_point <- "1h"
DMS_4h_mal_gostat_CC_hypo$time_point <- "4h"
DMS_24h_mal_gostat_CC_hypo$time_point <- "24h"
DMS_72h_mal_gostat_CC_hypo$time_point <- "72h"

DMS_mal_all_BP_hypo <- rbind(DMS_05h_mal_gostat_BP_hypo,DMS_1h_mal_gostat_BP_hypo,DMS_4h_mal_gostat_BP_hypo,DMS_24h_mal_gostat_BP_hypo,DMS_72h_mal_gostat_BP_hypo)
DMS_mal_all_MF_hypo <- rbind(DMS_05h_mal_gostat_MF_hypo,DMS_1h_mal_gostat_MF_hypo,DMS_4h_mal_gostat_MF_hypo,DMS_24h_mal_gostat_MF_hypo,DMS_72h_mal_gostat_MF_hypo)
DMS_mal_all_CC_hypo <- rbind(DMS_05h_mal_gostat_CC_hypo,DMS_1h_mal_gostat_CC_hypo,DMS_4h_mal_gostat_CC_hypo,DMS_24h_mal_gostat_CC_hypo,DMS_72h_mal_gostat_CC_hypo)

write.csv(DMS_mal_all_BP_hypo, "./shortTerm_exp/data/go_res_od/DMS_mal_combined_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMS_mal_all_MF_hypo, "./shortTerm_exp/data/go_res_od/DMS_mal_combined_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMS_mal_all_CC_hypo, "./shortTerm_exp/data/go_res_od/DMS_mal_combined_gostat_CC_hypo.csv", row.names = FALSE)

DMR_05h_fem_gostat_BP_hypo$time_point <- "0.5h"
DMR_1h_fem_gostat_BP_hypo$time_point <- "1h"
DMR_4h_fem_gostat_BP_hypo$time_point <- "4h"
DMR_24h_fem_gostat_BP_hypo$time_point <- "24h"
DMR_72h_fem_gostat_BP_hypo$time_point <- "72h"

DMR_05h_fem_gostat_MF_hypo$time_point <- "0.5h"
DMR_1h_fem_gostat_MF_hypo$time_point <- "1h"
DMR_4h_fem_gostat_MF_hypo$time_point <- "4h"
DMR_24h_fem_gostat_MF_hypo$time_point <- "24h"
DMR_72h_fem_gostat_MF_hypo$time_point <- "72h"

DMR_05h_fem_gostat_CC_hypo$time_point <- "0.5h"
DMR_1h_fem_gostat_CC_hypo$time_point <- "1h"
DMR_4h_fem_gostat_CC_hypo$time_point <- "4h"
DMR_24h_fem_gostat_CC_hypo$time_point <- "24h"
DMR_72h_fem_gostat_CC_hypo$time_point <- "72h"

DMR_fem_all_BP_hypo <- rbind(DMR_05h_fem_gostat_BP_hypo,DMR_1h_fem_gostat_BP_hypo,DMR_4h_fem_gostat_BP_hypo,DMR_24h_fem_gostat_BP_hypo,DMR_72h_fem_gostat_BP_hypo)
DMR_fem_all_MF_hypo <- rbind(DMR_05h_fem_gostat_MF_hypo,DMR_1h_fem_gostat_MF_hypo,DMR_4h_fem_gostat_MF_hypo,DMR_24h_fem_gostat_MF_hypo,DMR_72h_fem_gostat_MF_hypo)
DMR_fem_all_CC_hypo <- rbind(DMR_05h_fem_gostat_CC_hypo,DMR_1h_fem_gostat_CC_hypo,DMR_4h_fem_gostat_CC_hypo,DMR_24h_fem_gostat_CC_hypo,DMR_72h_fem_gostat_CC_hypo)

write.csv(DMR_fem_all_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_fem_combined_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_fem_all_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_fem_combined_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_fem_all_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_fem_combined_gostat_CC_hypo.csv", row.names = FALSE)

DMR_05h_mal_gostat_BP_hypo$time_point <- "0.5h"
DMR_1h_mal_gostat_BP_hypo$time_point <- "1h"
DMR_4h_mal_gostat_BP_hypo$time_point <- "4h"
DMR_24h_mal_gostat_BP_hypo$time_point <- "24h"
DMR_72h_mal_gostat_BP_hypo$time_point <- "72h"

DMR_05h_mal_gostat_MF_hypo$time_point <- "0.5h"
DMR_1h_mal_gostat_MF_hypo$time_point <- "1h"
DMR_4h_mal_gostat_MF_hypo$time_point <- "4h"
DMR_24h_mal_gostat_MF_hypo$time_point <- "24h"
DMR_72h_mal_gostat_MF_hypo$time_point <- "72h"

DMR_05h_mal_gostat_CC_hypo$time_point <- "0.5h"
DMR_1h_mal_gostat_CC_hypo$time_point <- "1h"
DMR_4h_mal_gostat_CC_hypo$time_point <- "4h"
DMR_24h_mal_gostat_CC_hypo$time_point <- "24h"
DMR_72h_mal_gostat_CC_hypo$time_point <- "72h"

DMR_mal_all_BP_hypo <- rbind(DMR_05h_mal_gostat_BP_hypo,DMR_1h_mal_gostat_BP_hypo,DMR_4h_mal_gostat_BP_hypo,DMR_24h_mal_gostat_BP_hypo,DMR_72h_mal_gostat_BP_hypo)
DMR_mal_all_MF_hypo <- rbind(DMR_05h_mal_gostat_MF_hypo,DMR_1h_mal_gostat_MF_hypo,DMR_4h_mal_gostat_MF_hypo,DMR_24h_mal_gostat_MF_hypo,DMR_72h_mal_gostat_MF_hypo)
DMR_mal_all_CC_hypo <- rbind(DMR_05h_mal_gostat_CC_hypo,DMR_1h_mal_gostat_CC_hypo,DMR_4h_mal_gostat_CC_hypo,DMR_24h_mal_gostat_CC_hypo,DMR_72h_mal_gostat_CC_hypo)

write.csv(DMR_mal_all_BP_hypo, "./shortTerm_exp/data/go_res_od/DMR_mal_combined_gostat_BP_hypo.csv", row.names = FALSE)
write.csv(DMR_mal_all_MF_hypo, "./shortTerm_exp/data/go_res_od/DMR_mal_combined_gostat_MF_hypo.csv", row.names = FALSE)
write.csv(DMR_mal_all_CC_hypo, "./shortTerm_exp/data/go_res_od/DMR_mal_combined_gostat_CC_hypo.csv", row.names = FALSE)
 
## plot go seq analysis ##
DMS.fem.BP <- plot.funct(DMS_05h_fem_gostat_BP, DMS_1h_fem_gostat_BP, DMS_4h_fem_gostat_BP, DMS_24h_fem_gostat_BP, DMS_72h_fem_gostat_BP,
               nrow(DMS_TSS_05h_fem), nrow(DMS_TSS_1h_fem), nrow(DMS_TSS_4h_fem), nrow(DMS_TSS_24h_fem), nrow(DMS_TSS_72h_fem),
               "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_fem_goplot_BP.tiff")

DMS.fem.MF <- plot.funct(DMS_05h_fem_gostat_MF, DMS_1h_fem_gostat_MF, DMS_4h_fem_gostat_MF, DMS_24h_fem_gostat_MF, DMS_72h_fem_gostat_MF,
           nrow(DMS_TSS_05h_fem), nrow(DMS_TSS_1h_fem), nrow(DMS_TSS_4h_fem), nrow(DMS_TSS_24h_fem), nrow(DMS_TSS_72h_fem),
           "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_fem_goplot_MF.tiff")

DMR.fem.BP <- plot.funct(DMR_05h_fem_gostat_BP, DMR_1h_fem_gostat_BP, DMR_4h_fem_gostat_BP, DMR_24h_fem_gostat_BP, DMR_72h_fem_gostat_BP,
           nrow(DMR_TSS_05h_fem), nrow(DMR_TSS_1h_fem), nrow(DMR_TSS_4h_fem), nrow(DMR_TSS_24h_fem), nrow(DMR_TSS_72h_fem),
           "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_fem_goplot_BP.tiff")

DMR.fem.MF <- plot.funct(DMR_05h_fem_gostat_MF, DMR_1h_fem_gostat_MF, DMR_4h_fem_gostat_MF, DMR_24h_fem_gostat_MF, DMR_72h_fem_gostat_MF,
           nrow(DMR_TSS_05h_fem), nrow(DMR_TSS_1h_fem), nrow(DMR_TSS_4h_fem), nrow(DMR_TSS_24h_fem), nrow(DMR_TSS_72h_fem),
           "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_fem_goplot_MF.tiff")

DMS.mal.BP <- plot.funct(DMS_05h_mal_gostat_BP, DMS_1h_mal_gostat_BP, DMS_4h_mal_gostat_BP, DMS_24h_mal_gostat_BP, DMS_72h_mal_gostat_BP,
           nrow(DMS_TSS_05h_mal), nrow(DMS_TSS_1h_mal), nrow(DMS_TSS_4h_mal), nrow(DMS_TSS_24h_mal), nrow(DMS_TSS_72h_mal),
           "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_mal_goplot_BP.tiff")

DMS.mal.MF <- plot.funct(DMS_05h_mal_gostat_MF, DMS_1h_mal_gostat_MF, DMS_4h_mal_gostat_MF, DMS_24h_mal_gostat_MF, DMS_72h_mal_gostat_MF,
           nrow(DMS_TSS_05h_mal), nrow(DMS_TSS_1h_mal), nrow(DMS_TSS_4h_mal), nrow(DMS_TSS_24h_mal), nrow(DMS_TSS_72h_mal),
           "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_mal_goplot_MF.tiff")

DMR.mal.BP <- plot.funct(DMR_05h_mal_gostat_BP, DMR_1h_mal_gostat_BP, DMR_4h_mal_gostat_BP, DMR_24h_mal_gostat_BP, DMR_72h_mal_gostat_BP,
           nrow(DMR_TSS_05h_mal), nrow(DMR_TSS_1h_mal), nrow(DMR_TSS_4h_mal), nrow(DMR_TSS_24h_mal), nrow(DMR_TSS_72h_mal),
           "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_mal_goplot_BP.tiff")

DMR.mal.MF <- plot.funct(DMR_05h_mal_gostat_MF, DMR_1h_mal_gostat_MF, DMR_4h_mal_gostat_MF, DMR_24h_mal_gostat_MF, DMR_72h_mal_gostat_MF,
           nrow(DMR_TSS_05h_mal), nrow(DMR_TSS_1h_mal), nrow(DMR_TSS_4h_mal), nrow(DMR_TSS_24h_mal), nrow(DMR_TSS_72h_mal),
           "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_mal_goplot_MF.tiff")

#hyper
DMS.fem.BP <- plot.funct(DMS_05h_fem_gostat_BP_hyper, DMS_1h_fem_gostat_BP_hyper, DMS_4h_fem_gostat_BP_hyper, DMS_24h_fem_gostat_BP_hyper, DMS_72h_fem_gostat_BP_hyper,
                         nrow(DMS_TSS_05h_fem_hyper), nrow(DMS_TSS_1h_fem_hyper), nrow(DMS_TSS_4h_fem_hyper), nrow(DMS_TSS_24h_fem_hyper), nrow(DMS_TSS_72h_fem_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_fem_goplot_BP_hyper.tiff")

DMS.fem.MF <- plot.funct(DMS_05h_fem_gostat_MF_hyper, DMS_1h_fem_gostat_MF_hyper, DMS_4h_fem_gostat_MF_hyper, DMS_24h_fem_gostat_MF_hyper, DMS_72h_fem_gostat_MF_hyper,
                         nrow(DMS_TSS_05h_fem_hyper), nrow(DMS_TSS_1h_fem_hyper), nrow(DMS_TSS_4h_fem_hyper), nrow(DMS_TSS_24h_fem_hyper), nrow(DMS_TSS_72h_fem_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_fem_goplot_MF_hyper.tiff")

DMR.fem.BP <- plot.funct(DMR_05h_fem_gostat_BP_hyper, DMR_1h_fem_gostat_BP_hyper, DMR_4h_fem_gostat_BP_hyper, DMR_24h_fem_gostat_BP_hyper, DMR_72h_fem_gostat_BP_hyper,
                         nrow(DMR_TSS_05h_fem_hyper), nrow(DMR_TSS_1h_fem_hyper), nrow(DMR_TSS_4h_fem_hyper), nrow(DMR_TSS_24h_fem_hyper), nrow(DMR_TSS_72h_fem_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_fem_goplot_BP_hyper.tiff")

DMR.fem.MF <- plot.funct(DMR_05h_fem_gostat_MF_hyper, DMR_1h_fem_gostat_MF_hyper, DMR_4h_fem_gostat_MF_hyper, DMR_24h_fem_gostat_MF_hyper, DMR_72h_fem_gostat_MF_hyper,
                         nrow(DMR_TSS_05h_fem_hyper), nrow(DMR_TSS_1h_fem_hyper), nrow(DMR_TSS_4h_fem_hyper), nrow(DMR_TSS_24h_fem_hyper), nrow(DMR_TSS_72h_fem_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_fem_goplot_MF_hyper.tiff")

DMS.mal.BP <- plot.funct(DMS_05h_mal_gostat_BP_hyper, DMS_1h_mal_gostat_BP_hyper, DMS_4h_mal_gostat_BP_hyper, DMS_24h_mal_gostat_BP_hyper, DMS_72h_mal_gostat_BP_hyper,
                         nrow(DMS_TSS_05h_mal_hyper), nrow(DMS_TSS_1h_mal_hyper), nrow(DMS_TSS_4h_mal_hyper), nrow(DMS_TSS_24h_mal_hyper), nrow(DMS_TSS_72h_mal_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_mal_goplot_BP_hyper.tiff")

DMS.mal.MF <- plot.funct(DMS_05h_mal_gostat_MF_hyper, DMS_1h_mal_gostat_MF_hyper, DMS_4h_mal_gostat_MF_hyper, DMS_24h_mal_gostat_MF_hyper, DMS_72h_mal_gostat_MF_hyper,
                         nrow(DMS_TSS_05h_mal_hyper), nrow(DMS_TSS_1h_mal_hyper), nrow(DMS_TSS_4h_mal_hyper), nrow(DMS_TSS_24h_mal_hyper), nrow(DMS_TSS_72h_mal_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_mal_goplot_MF_hyper.tiff")

DMR.mal.BP <- plot.funct(DMR_05h_mal_gostat_BP_hyper, DMR_1h_mal_gostat_BP_hyper, DMR_4h_mal_gostat_BP_hyper, DMR_24h_mal_gostat_BP_hyper, DMR_72h_mal_gostat_BP_hyper,
                         nrow(DMR_TSS_05h_mal_hyper), nrow(DMR_TSS_1h_mal_hyper), nrow(DMR_TSS_4h_mal_hyper), nrow(DMR_TSS_24h_mal_hyper), nrow(DMR_TSS_72h_mal_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_mal_goplot_BP_hyper.tiff")

DMR.mal.MF <- plot.funct(DMR_05h_mal_gostat_MF_hyper, DMR_1h_mal_gostat_MF_hyper, DMR_4h_mal_gostat_MF_hyper, DMR_24h_mal_gostat_MF_hyper, DMR_72h_mal_gostat_MF_hyper,
                         nrow(DMR_TSS_05h_mal_hyper), nrow(DMR_TSS_1h_mal_hyper), nrow(DMR_TSS_4h_mal_hyper), nrow(DMR_TSS_24h_mal_hyper), nrow(DMR_TSS_72h_mal_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_mal_goplot_MF_hyper.tiff")

#hypo
DMS.fem.BP <- plot.funct(DMS_05h_fem_gostat_BP_hypo, DMS_1h_fem_gostat_BP_hypo, DMS_4h_fem_gostat_BP_hypo, DMS_24h_fem_gostat_BP_hypo, DMS_72h_fem_gostat_BP_hypo,
                         nrow(DMS_TSS_05h_fem_hypo), nrow(DMS_TSS_1h_fem_hypo), nrow(DMS_TSS_4h_fem_hypo), nrow(DMS_TSS_24h_fem_hypo), nrow(DMS_TSS_72h_fem_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_fem_goplot_BP_hypo.tiff")

DMS.fem.MF <- plot.funct(DMS_05h_fem_gostat_MF_hypo, DMS_1h_fem_gostat_MF_hypo, DMS_4h_fem_gostat_MF_hypo, DMS_24h_fem_gostat_MF_hypo, DMS_72h_fem_gostat_MF_hypo,
                         nrow(DMS_TSS_05h_fem_hypo), nrow(DMS_TSS_1h_fem_hypo), nrow(DMS_TSS_4h_fem_hypo), nrow(DMS_TSS_24h_fem_hypo), nrow(DMS_TSS_72h_fem_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_fem_goplot_MF_hypo.tiff")

DMR.fem.BP <- plot.funct(DMR_05h_fem_gostat_BP_hypo, DMR_1h_fem_gostat_BP_hypo, DMR_4h_fem_gostat_BP_hypo, DMR_24h_fem_gostat_BP_hypo, DMR_72h_fem_gostat_BP_hypo,
                         nrow(DMR_TSS_05h_fem_hypo), nrow(DMR_TSS_1h_fem_hypo), nrow(DMR_TSS_4h_fem_hypo), nrow(DMR_TSS_24h_fem_hypo), nrow(DMR_TSS_72h_fem_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_fem_goplot_BP_hypo.tiff")

DMR.fem.MF <- plot.funct(DMR_05h_fem_gostat_MF_hypo, DMR_1h_fem_gostat_MF_hypo, DMR_4h_fem_gostat_MF_hypo, DMR_24h_fem_gostat_MF_hypo, DMR_72h_fem_gostat_MF_hypo,
                         nrow(DMR_TSS_05h_fem_hypo), nrow(DMR_TSS_1h_fem_hypo), nrow(DMR_TSS_4h_fem_hypo), nrow(DMR_TSS_24h_fem_hypo), nrow(DMR_TSS_72h_fem_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_fem_goplot_MF_hypo.tiff")

DMS.mal.BP <- plot.funct(DMS_05h_mal_gostat_BP_hypo, DMS_1h_mal_gostat_BP_hypo, DMS_4h_mal_gostat_BP_hypo, DMS_24h_mal_gostat_BP_hypo, DMS_72h_mal_gostat_BP_hypo,
                         nrow(DMS_TSS_05h_mal_hypo), nrow(DMS_TSS_1h_mal_hypo), nrow(DMS_TSS_4h_mal_hypo), nrow(DMS_TSS_24h_mal_hypo), nrow(DMS_TSS_72h_mal_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_mal_goplot_BP_hypo.tiff")

DMS.mal.MF <- plot.funct(DMS_05h_mal_gostat_MF_hypo, DMS_1h_mal_gostat_MF_hypo, DMS_4h_mal_gostat_MF_hypo, DMS_24h_mal_gostat_MF_hypo, DMS_72h_mal_gostat_MF_hypo,
                         nrow(DMS_TSS_05h_mal_hypo), nrow(DMS_TSS_1h_mal_hypo), nrow(DMS_TSS_4h_mal_hypo), nrow(DMS_TSS_24h_mal_hypo), nrow(DMS_TSS_72h_mal_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_mal_goplot_MF_hypo.tiff")

DMR.mal.BP <- plot.funct(DMR_05h_mal_gostat_BP_hypo, DMR_1h_mal_gostat_BP_hypo, DMR_4h_mal_gostat_BP_hypo, DMR_24h_mal_gostat_BP_hypo, DMR_72h_mal_gostat_BP_hypo,
                         nrow(DMR_TSS_05h_mal_hypo), nrow(DMR_TSS_1h_mal_hypo), nrow(DMR_TSS_4h_mal_hypo), nrow(DMR_TSS_24h_mal_hypo), nrow(DMR_TSS_72h_mal_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_mal_goplot_BP_hypo.tiff")

DMR.mal.MF <- plot.funct(DMR_05h_mal_gostat_MF_hypo, DMR_1h_mal_gostat_MF_hypo, DMR_4h_mal_gostat_MF_hypo, DMR_24h_mal_gostat_MF_hypo, DMR_72h_mal_gostat_MF_hypo,
                         nrow(DMR_TSS_05h_mal_hypo), nrow(DMR_TSS_1h_mal_hypo), nrow(DMR_TSS_4h_mal_hypo), nrow(DMR_TSS_24h_mal_hypo), nrow(DMR_TSS_72h_mal_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_mal_goplot_MF_hypo.tiff")

#plot go seq analysis
DMS.fem.BP <- plot.funct.top5(DMS_05h_fem_gostat_BP, DMS_1h_fem_gostat_BP, DMS_4h_fem_gostat_BP, DMS_24h_fem_gostat_BP, DMS_72h_fem_gostat_BP,
                         nrow(DMS_TSS_05h_fem), nrow(DMS_TSS_1h_fem), nrow(DMS_TSS_4h_fem), nrow(DMS_TSS_24h_fem), nrow(DMS_TSS_72h_fem),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_fem_goplot_BP.tiff")

DMS.fem.MF <- plot.funct.top5(DMS_05h_fem_gostat_MF, DMS_1h_fem_gostat_MF, DMS_4h_fem_gostat_MF, DMS_24h_fem_gostat_MF, DMS_72h_fem_gostat_MF,
                         nrow(DMS_TSS_05h_fem), nrow(DMS_TSS_1h_fem), nrow(DMS_TSS_4h_fem), nrow(DMS_TSS_24h_fem), nrow(DMS_TSS_72h_fem),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_fem_goplot_MF.tiff")

DMR.fem.BP <- plot.funct.top5(DMR_05h_fem_gostat_BP, DMR_1h_fem_gostat_BP, DMR_4h_fem_gostat_BP, DMR_24h_fem_gostat_BP, DMR_72h_fem_gostat_BP,
                         nrow(DMR_TSS_05h_fem), nrow(DMR_TSS_1h_fem), nrow(DMR_TSS_4h_fem), nrow(DMR_TSS_24h_fem), nrow(DMR_TSS_72h_fem),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_fem_goplot_BP.tiff")

DMR.fem.MF <- plot.funct.top5(DMR_05h_fem_gostat_MF, DMR_1h_fem_gostat_MF, DMR_4h_fem_gostat_MF, DMR_24h_fem_gostat_MF, DMR_72h_fem_gostat_MF,
                         nrow(DMR_TSS_05h_fem), nrow(DMR_TSS_1h_fem), nrow(DMR_TSS_4h_fem), nrow(DMR_TSS_24h_fem), nrow(DMR_TSS_72h_fem),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_fem_goplot_MF.tiff")

DMS.mal.BP <- plot.funct.top5(DMS_05h_mal_gostat_BP, DMS_1h_mal_gostat_BP, DMS_4h_mal_gostat_BP, DMS_24h_mal_gostat_BP, DMS_72h_mal_gostat_BP,
                         nrow(DMS_TSS_05h_mal), nrow(DMS_TSS_1h_mal), nrow(DMS_TSS_4h_mal), nrow(DMS_TSS_24h_mal), nrow(DMS_TSS_72h_mal),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_mal_goplot_BP.tiff")

DMS.mal.MF <- plot.funct.top5(DMS_05h_mal_gostat_MF, DMS_1h_mal_gostat_MF, DMS_4h_mal_gostat_MF, DMS_24h_mal_gostat_MF, DMS_72h_mal_gostat_MF,
                         nrow(DMS_TSS_05h_mal), nrow(DMS_TSS_1h_mal), nrow(DMS_TSS_4h_mal), nrow(DMS_TSS_24h_mal), nrow(DMS_TSS_72h_mal),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_mal_goplot_MF.tiff")

DMR.mal.BP <- plot.funct.top5(DMR_05h_mal_gostat_BP, DMR_1h_mal_gostat_BP, DMR_4h_mal_gostat_BP, DMR_24h_mal_gostat_BP, DMR_72h_mal_gostat_BP,
                         nrow(DMR_TSS_05h_mal), nrow(DMR_TSS_1h_mal), nrow(DMR_TSS_4h_mal), nrow(DMR_TSS_24h_mal), nrow(DMR_TSS_72h_mal),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_mal_goplot_BP.tiff")

DMR.mal.MF <- plot.funct.top5(DMR_05h_mal_gostat_MF, DMR_1h_mal_gostat_MF, DMR_4h_mal_gostat_MF, DMR_24h_mal_gostat_MF, DMR_72h_mal_gostat_MF,
                         nrow(DMR_TSS_05h_mal), nrow(DMR_TSS_1h_mal), nrow(DMR_TSS_4h_mal), nrow(DMR_TSS_24h_mal), nrow(DMR_TSS_72h_mal),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_mal_goplot_MF.tiff")

#hyper
DMS.fem.BP.hyper <- plot.funct.top5(DMS_05h_fem_gostat_BP_hyper, DMS_1h_fem_gostat_BP_hyper, DMS_4h_fem_gostat_BP_hyper, DMS_24h_fem_gostat_BP_hyper, DMS_72h_fem_gostat_BP_hyper,
                         nrow(DMS_TSS_05h_fem_hyper), nrow(DMS_TSS_1h_fem_hyper), nrow(DMS_TSS_4h_fem_hyper), nrow(DMS_TSS_24h_fem_hyper), nrow(DMS_TSS_72h_fem_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_fem_goplot_BP_hyper.tiff")

DMS.fem.MF.hyper <- plot.funct.top5(DMS_05h_fem_gostat_MF_hyper, DMS_1h_fem_gostat_MF_hyper, DMS_4h_fem_gostat_MF_hyper, DMS_24h_fem_gostat_MF_hyper, DMS_72h_fem_gostat_MF_hyper,
                         nrow(DMS_TSS_05h_fem_hyper), nrow(DMS_TSS_1h_fem_hyper), nrow(DMS_TSS_4h_fem_hyper), nrow(DMS_TSS_24h_fem_hyper), nrow(DMS_TSS_72h_fem_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_fem_goplot_MF_hyper.tiff")

DMR.fem.BP.hyper <- plot.funct.top5(DMR_05h_fem_gostat_BP_hyper, DMR_1h_fem_gostat_BP_hyper, DMR_4h_fem_gostat_BP_hyper, DMR_24h_fem_gostat_BP_hyper, DMR_72h_fem_gostat_BP_hyper,
                         nrow(DMR_TSS_05h_fem_hyper), nrow(DMR_TSS_1h_fem_hyper), nrow(DMR_TSS_4h_fem_hyper), nrow(DMR_TSS_24h_fem_hyper), nrow(DMR_TSS_72h_fem_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_fem_goplot_BP_hyper.tiff")

DMR.fem.MF.hyper <- plot.funct.top5(DMR_05h_fem_gostat_MF_hyper, DMR_1h_fem_gostat_MF_hyper, DMR_4h_fem_gostat_MF_hyper, DMR_24h_fem_gostat_MF_hyper, DMR_72h_fem_gostat_MF_hyper,
                         nrow(DMR_TSS_05h_fem_hyper), nrow(DMR_TSS_1h_fem_hyper), nrow(DMR_TSS_4h_fem_hyper), nrow(DMR_TSS_24h_fem_hyper), nrow(DMR_TSS_72h_fem_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_fem_goplot_MF_hyper.tiff")

DMS.mal.BP.hyper <- plot.funct.top5(DMS_05h_mal_gostat_BP_hyper, DMS_1h_mal_gostat_BP_hyper, DMS_4h_mal_gostat_BP_hyper, DMS_24h_mal_gostat_BP_hyper, DMS_72h_mal_gostat_BP_hyper,
                         nrow(DMS_TSS_05h_mal_hyper), nrow(DMS_TSS_1h_mal_hyper), nrow(DMS_TSS_4h_mal_hyper), nrow(DMS_TSS_24h_mal_hyper), nrow(DMS_TSS_72h_mal_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_mal_goplot_BP_hyper.tiff")

DMS.mal.MF.hyper <- plot.funct.top5(DMS_05h_mal_gostat_MF_hyper, DMS_1h_mal_gostat_MF_hyper, DMS_4h_mal_gostat_MF_hyper, DMS_24h_mal_gostat_MF_hyper, DMS_72h_mal_gostat_MF_hyper,
                         nrow(DMS_TSS_05h_mal_hyper), nrow(DMS_TSS_1h_mal_hyper), nrow(DMS_TSS_4h_mal_hyper), nrow(DMS_TSS_24h_mal_hyper), nrow(DMS_TSS_72h_mal_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_mal_goplot_MF_hyper.tiff")

DMR.mal.BP.hyper <- plot.funct.top5(DMR_05h_mal_gostat_BP_hyper, DMR_1h_mal_gostat_BP_hyper, DMR_4h_mal_gostat_BP_hyper, DMR_24h_mal_gostat_BP_hyper, DMR_72h_mal_gostat_BP_hyper,
                         nrow(DMR_TSS_05h_mal_hyper), nrow(DMR_TSS_1h_mal_hyper), nrow(DMR_TSS_4h_mal_hyper), nrow(DMR_TSS_24h_mal_hyper), nrow(DMR_TSS_72h_mal_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_mal_goplot_BP_hyper.tiff")

DMR.mal.MF.hyper <- plot.funct.top5(DMR_05h_mal_gostat_MF_hyper, DMR_1h_mal_gostat_MF_hyper, DMR_4h_mal_gostat_MF_hyper, DMR_24h_mal_gostat_MF_hyper, DMR_72h_mal_gostat_MF_hyper,
                         nrow(DMR_TSS_05h_mal_hyper), nrow(DMR_TSS_1h_mal_hyper), nrow(DMR_TSS_4h_mal_hyper), nrow(DMR_TSS_24h_mal_hyper), nrow(DMR_TSS_72h_mal_hyper),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_mal_goplot_MF_hyper.tiff")

#hypo
DMS.fem.BP.hypo <- plot.funct.top5(DMS_05h_fem_gostat_BP_hypo, DMS_1h_fem_gostat_BP_hypo, DMS_4h_fem_gostat_BP_hypo, DMS_24h_fem_gostat_BP_hypo, DMS_72h_fem_gostat_BP_hypo,
                         nrow(DMS_TSS_05h_fem_hypo), nrow(DMS_TSS_1h_fem_hypo), nrow(DMS_TSS_4h_fem_hypo), nrow(DMS_TSS_24h_fem_hypo), nrow(DMS_TSS_72h_fem_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_fem_goplot_BP_hypo.tiff")

DMS.fem.MF.hypo <- plot.funct.top5(DMS_05h_fem_gostat_MF_hypo, DMS_1h_fem_gostat_MF_hypo, DMS_4h_fem_gostat_MF_hypo, DMS_24h_fem_gostat_MF_hypo, DMS_72h_fem_gostat_MF_hypo,
                         nrow(DMS_TSS_05h_fem_hypo), nrow(DMS_TSS_1h_fem_hypo), nrow(DMS_TSS_4h_fem_hypo), nrow(DMS_TSS_24h_fem_hypo), nrow(DMS_TSS_72h_fem_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_fem_goplot_MF_hypo.tiff")

DMR.fem.BP.hypo <- plot.funct.top5(DMR_05h_fem_gostat_BP_hypo, DMR_1h_fem_gostat_BP_hypo, DMR_4h_fem_gostat_BP_hypo, DMR_24h_fem_gostat_BP_hypo, DMR_72h_fem_gostat_BP_hypo,
                         nrow(DMR_TSS_05h_fem_hypo), nrow(DMR_TSS_1h_fem_hypo), nrow(DMR_TSS_4h_fem_hypo), nrow(DMR_TSS_24h_fem_hypo), nrow(DMR_TSS_72h_fem_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_fem_goplot_BP_hypo.tiff")

DMR.fem.MF.hypo <- plot.funct.top5(DMR_05h_fem_gostat_MF_hypo, DMR_1h_fem_gostat_MF_hypo, DMR_4h_fem_gostat_MF_hypo, DMR_24h_fem_gostat_MF_hypo, DMR_72h_fem_gostat_MF_hypo,
                         nrow(DMR_TSS_05h_fem_hypo), nrow(DMR_TSS_1h_fem_hypo), nrow(DMR_TSS_4h_fem_hypo), nrow(DMR_TSS_24h_fem_hypo), nrow(DMR_TSS_72h_fem_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_fem_goplot_MF_hypo.tiff")

DMS.mal.BP.hypo <- plot.funct.top5(DMS_05h_mal_gostat_BP_hypo, DMS_1h_mal_gostat_BP_hypo, DMS_4h_mal_gostat_BP_hypo, DMS_24h_mal_gostat_BP_hypo, DMS_72h_mal_gostat_BP_hypo,
                         nrow(DMS_TSS_05h_mal_hypo), nrow(DMS_TSS_1h_mal_hypo), nrow(DMS_TSS_4h_mal_hypo), nrow(DMS_TSS_24h_mal_hypo), nrow(DMS_TSS_72h_mal_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_mal_goplot_BP_hypo.tiff")

DMS.mal.MF.hypo <- plot.funct.top5(DMS_05h_mal_gostat_MF_hypo, DMS_1h_mal_gostat_MF_hypo, DMS_4h_mal_gostat_MF_hypo, DMS_24h_mal_gostat_MF_hypo, DMS_72h_mal_gostat_MF_hypo,
                         nrow(DMS_TSS_05h_mal_hypo), nrow(DMS_TSS_1h_mal_hypo), nrow(DMS_TSS_4h_mal_hypo), nrow(DMS_TSS_24h_mal_hypo), nrow(DMS_TSS_72h_mal_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_mal_goplot_MF_hypo.tiff")

DMR.mal.BP.hypo <- plot.funct.top5(DMR_05h_mal_gostat_BP_hypo, DMR_1h_mal_gostat_BP_hypo, DMR_4h_mal_gostat_BP_hypo, DMR_24h_mal_gostat_BP_hypo, DMR_72h_mal_gostat_BP_hypo,
                         nrow(DMR_TSS_05h_mal_hypo), nrow(DMR_TSS_1h_mal_hypo), nrow(DMR_TSS_4h_mal_hypo), nrow(DMR_TSS_24h_mal_hypo), nrow(DMR_TSS_72h_mal_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_mal_goplot_BP_hypo.tiff")

DMR.mal.MF.hypo <- plot.funct.top5(DMR_05h_mal_gostat_MF_hypo, DMR_1h_mal_gostat_MF_hypo, DMR_4h_mal_gostat_MF_hypo, DMR_24h_mal_gostat_MF_hypo, DMR_72h_mal_gostat_MF_hypo,
                         nrow(DMR_TSS_05h_mal_hypo), nrow(DMR_TSS_1h_mal_hypo), nrow(DMR_TSS_4h_mal_hypo), nrow(DMR_TSS_24h_mal_hypo), nrow(DMR_TSS_72h_mal_hypo),
                         "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_mal_goplot_MF_hypo.tiff")

#make panels
DMR.fem.panel <- ggarrange(DMR.fem.BP.hyper, DMR.fem.BP.hypo, nrow = 2, ncol = 1, labels = c("A", "B"), common.legend = TRUE, legend = "right")
ggsave(plot = DMR.fem.panel, filename = "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_fem_goplot_BP_panel.tiff", 
       width = 8, height = 11, units = "in", dpi = 300)

DMR.mal.panel <- ggarrange(DMR.mal.BP.hyper, DMR.mal.BP.hypo, nrow = 2, ncol = 1, labels = c("A", "B"), common.legend = TRUE, legend = "right")
ggsave(plot = DMR.mal.panel, filename = "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMR_mal_goplot_BP_panel.tiff", 
       width = 8, height = 11, units = "in", dpi = 300)

DMS.fem.panel <- ggarrange(DMS.fem.BP.hyper, DMS.fem.BP.hypo, nrow = 2, ncol = 1, labels = c("A", "B"), common.legend = TRUE, legend = "right")
ggsave(plot = DMS.fem.panel, filename = "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_fem_goplot_BP_panel.tiff", 
       width = 8, height = 11, units = "in", dpi = 300)

DMS.mal.panel <- ggarrange(DMS.mal.BP.hyper, DMS.mal.BP.hypo, nrow = 2, ncol = 1, labels = c("A", "B"), common.legend = TRUE, legend = "right")
ggsave(plot = DMS.mal.panel, filename = "./shortTerm_exp/plots/finalized_tiff/go_plots_od/top5/DMS_mal_goplot_BP_panel.tiff", 
       width = 8, height = 11, units = "in", dpi = 300)

# plot results of pooled analysis 
plot_pooled(DMR_all_fem_gostat_BP_hyper, nrow(DMR_TSS_all_fem_hyper),
            DMR_all_fem_gostat_BP_hypo, nrow(DMR_TSS_all_fem_hypo),
            "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_fem_pooled_goplot_BP.tiff")

plot_pooled(DMS_all_fem_gostat_BP_hyper, nrow(DMS_TSS_all_fem_hyper),
            DMS_all_fem_gostat_BP_hypo, nrow(DMS_TSS_all_fem_hypo),
            "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_fem_pooled_goplot_BP.tiff")

plot_pooled(DMR_all_mal_gostat_BP_hyper, nrow(DMR_TSS_all_mal_hyper),
            DMR_all_mal_gostat_BP_hypo, nrow(DMR_TSS_all_mal_hypo),
            "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMR_mal_pooled_goplot_BP.tiff")

plot_pooled(DMS_all_mal_gostat_BP_hyper, nrow(DMS_TSS_all_mal_hyper),
            DMS_all_mal_gostat_BP_hypo, nrow(DMS_TSS_all_mal_hypo),
            "./shortTerm_exp/plots/finalized_tiff/go_plots_od/DMS_mal_pooled_goplot_BP.tiff")

## Get lists of differentially methylated genes ##
gene.list.funct(DMS_TSS_all_fem, "./shortTerm_exp/data/gene_lists_od/DMS_allTP_fem_genelist.csv")
gene.list.funct(DMS_TSS_05h_fem, "./shortTerm_exp/data/gene_lists_od/DMS_05h_fem_genelist.csv")
gene.list.funct(DMS_TSS_1h_fem, "./shortTerm_exp/data/gene_lists_od/DMS_1h_fem_genelist.csv")
gene.list.funct(DMS_TSS_4h_fem, "./shortTerm_exp/data/gene_lists_od/DMS_4h_fem_genelist.csv")
gene.list.funct(DMS_TSS_24h_fem, "./shortTerm_exp/data/gene_lists_od/DMS_24h_fem_genelist.csv")
gene.list.funct(DMS_TSS_72h_fem, "./shortTerm_exp/data/gene_lists_od/DMS_72h_fem_genelist.csv")

gene.list.funct(DMS_TSS_all_mal, "./shortTerm_exp/data/gene_lists_od/DMS_allTP_mal_genelist.csv")
gene.list.funct(DMS_TSS_05h_mal, "./shortTerm_exp/data/gene_lists_od/DMS_05h_mal_genelist.csv")
gene.list.funct(DMS_TSS_1h_mal, "./shortTerm_exp/data/gene_lists_od/DMS_1h_mal_genelist.csv")
gene.list.funct(DMS_TSS_4h_mal, "./shortTerm_exp/data/gene_lists_od/DMS_4h_mal_genelist.csv")
gene.list.funct(DMS_TSS_24h_mal, "./shortTerm_exp/data/gene_lists_od/DMS_24h_mal_genelist.csv")
gene.list.funct(DMS_TSS_72h_mal, "./shortTerm_exp/data/gene_lists_od/DMS_72h_mal_genelist.csv")

gene.list.funct(DMR_TSS_all_fem, "./shortTerm_exp/data/gene_lists_od/DMR_allTP_fem_genelist.csv")
gene.list.funct(DMR_TSS_05h_fem, "./shortTerm_exp/data/gene_lists_od/DMR_05h_fem_genelist.csv")
gene.list.funct(DMR_TSS_1h_fem, "./shortTerm_exp/data/gene_lists_od/DMR_1h_fem_genelist.csv")
gene.list.funct(DMR_TSS_4h_fem, "./shortTerm_exp/data/gene_lists_od/DMR_4h_fem_genelist.csv")
gene.list.funct(DMR_TSS_24h_fem, "./shortTerm_exp/data/gene_lists_od/DMR_24h_fem_genelist.csv")
gene.list.funct(DMR_TSS_72h_fem, "./shortTerm_exp/data/gene_lists_od/DMR_72h_fem_genelist.csv")

gene.list.funct(DMR_TSS_all_mal, "./shortTerm_exp/data/gene_lists_od/DMR_allTP_mal_genelist.csv")
gene.list.funct(DMR_TSS_05h_mal, "./shortTerm_exp/data/gene_lists_od/DMR_05h_mal_genelist.csv")
gene.list.funct(DMR_TSS_1h_mal, "./shortTerm_exp/data/gene_lists_od/DMR_1h_mal_genelist.csv")
gene.list.funct(DMR_TSS_4h_mal, "./shortTerm_exp/data/gene_lists_od/DMR_4h_mal_genelist.csv")
gene.list.funct(DMR_TSS_24h_mal, "./shortTerm_exp/data/gene_lists_od/DMR_24h_mal_genelist.csv")
gene.list.funct(DMR_TSS_72h_mal, "./shortTerm_exp/data/gene_lists_od/DMR_72h_mal_genelist.csv")

#hyper
gene.list.funct(DMS_TSS_all_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_allTP_fem_genelist.csv")
gene.list.funct(DMS_TSS_05h_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_05h_fem_genelist.csv")
gene.list.funct(DMS_TSS_1h_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_1h_fem_genelist.csv")
gene.list.funct(DMS_TSS_4h_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_4h_fem_genelist.csv")
gene.list.funct(DMS_TSS_24h_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_24h_fem_genelist.csv")
gene.list.funct(DMS_TSS_72h_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_72h_fem_genelist.csv")

gene.list.funct(DMS_TSS_all_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_allTP_mal_genelist.csv")
gene.list.funct(DMS_TSS_05h_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_05h_mal_genelist.csv")
gene.list.funct(DMS_TSS_1h_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_1h_mal_genelist.csv")
gene.list.funct(DMS_TSS_4h_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_4h_mal_genelist.csv")
gene.list.funct(DMS_TSS_24h_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_24h_mal_genelist.csv")
gene.list.funct(DMS_TSS_72h_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMS_72h_mal_genelist.csv")

gene.list.funct(DMR_TSS_all_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_allTP_fem_genelist.csv")
gene.list.funct(DMR_TSS_05h_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_05h_fem_genelist.csv")
gene.list.funct(DMR_TSS_1h_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_1h_fem_genelist.csv")
gene.list.funct(DMR_TSS_4h_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_4h_fem_genelist.csv")
gene.list.funct(DMR_TSS_24h_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_24h_fem_genelist.csv")
gene.list.funct(DMR_TSS_72h_fem_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_72h_fem_genelist.csv")

gene.list.funct(DMR_TSS_all_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_allTP_mal_genelist.csv")
gene.list.funct(DMR_TSS_05h_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_05h_mal_genelist.csv")
gene.list.funct(DMR_TSS_1h_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_1h_mal_genelist.csv")
gene.list.funct(DMR_TSS_4h_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_4h_mal_genelist.csv")
gene.list.funct(DMR_TSS_24h_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_24h_mal_genelist.csv")
gene.list.funct(DMR_TSS_72h_mal_hyper,  "./shortTerm_exp/data/gene_lists_od/DMR_72h_mal_genelist.csv")

#hypo
gene.list.funct(DMS_TSS_all_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_allTP_fem_genelist.csv")
gene.list.funct(DMS_TSS_05h_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_05h_fem_genelist.csv")
gene.list.funct(DMS_TSS_1h_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_1h_fem_genelist.csv")
gene.list.funct(DMS_TSS_4h_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_4h_fem_genelist.csv")
gene.list.funct(DMS_TSS_24h_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_24h_fem_genelist.csv")
gene.list.funct(DMS_TSS_72h_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_72h_fem_genelist.csv")

gene.list.funct(DMS_TSS_all_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_allTP_mal_genelist.csv")
gene.list.funct(DMS_TSS_05h_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_05h_mal_genelist.csv")
gene.list.funct(DMS_TSS_1h_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_1h_mal_genelist.csv")
gene.list.funct(DMS_TSS_4h_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_4h_mal_genelist.csv")
gene.list.funct(DMS_TSS_24h_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_24h_mal_genelist.csv")
gene.list.funct(DMS_TSS_72h_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMS_72h_mal_genelist.csv")

gene.list.funct(DMR_TSS_all_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_allTP_fem_genelist.csv")
gene.list.funct(DMR_TSS_05h_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_05h_fem_genelist.csv")
gene.list.funct(DMR_TSS_1h_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_1h_fem_genelist.csv")
gene.list.funct(DMR_TSS_4h_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_4h_fem_genelist.csv")
gene.list.funct(DMR_TSS_24h_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_24h_fem_genelist.csv")
gene.list.funct(DMR_TSS_72h_fem_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_72h_fem_genelist.csv")

gene.list.funct(DMR_TSS_all_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_allTP_mal_genelist.csv")
gene.list.funct(DMR_TSS_05h_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_05h_mal_genelist.csv")
gene.list.funct(DMR_TSS_1h_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_1h_mal_genelist.csv")
gene.list.funct(DMR_TSS_4h_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_4h_mal_genelist.csv")
gene.list.funct(DMR_TSS_24h_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_24h_mal_genelist.csv")
gene.list.funct(DMR_TSS_72h_mal_hypo,  "./shortTerm_exp/data/gene_lists_od/DMR_72h_mal_genelist.csv")

#check ranges
max(nrow(DMR_05h_fem_gostat_BP_hyper), nrow(DMR_1h_fem_gostat_BP_hyper), nrow(DMR_4h_fem_gostat_BP_hyper),
    nrow(DMR_24h_fem_gostat_BP_hyper), nrow(DMR_72h_fem_gostat_BP_hyper), 
    nrow(DMR_05h_fem_gostat_BP_hypo), nrow(DMR_1h_fem_gostat_BP_hypo), nrow(DMR_4h_fem_gostat_BP_hypo),
    nrow(DMR_24h_fem_gostat_BP_hypo), nrow(DMR_72h_fem_gostat_BP_hypo))

min(nrow(DMR_05h_fem_gostat_BP_hyper), nrow(DMR_1h_fem_gostat_BP_hyper), nrow(DMR_4h_fem_gostat_BP_hyper),
    nrow(DMR_24h_fem_gostat_BP_hyper), nrow(DMR_72h_fem_gostat_BP_hyper), 
    nrow(DMR_05h_fem_gostat_BP_hypo), nrow(DMR_1h_fem_gostat_BP_hypo), nrow(DMR_4h_fem_gostat_BP_hypo),
    nrow(DMR_24h_fem_gostat_BP_hypo), nrow(DMR_72h_fem_gostat_BP_hypo))

max(nrow(DMS_05h_fem_gostat_BP_hyper), nrow(DMS_1h_fem_gostat_BP_hyper), nrow(DMS_4h_fem_gostat_BP_hyper),
    nrow(DMS_24h_fem_gostat_BP_hyper), nrow(DMS_72h_fem_gostat_BP_hyper), 
    nrow(DMS_05h_fem_gostat_BP_hypo), nrow(DMS_1h_fem_gostat_BP_hypo), nrow(DMS_4h_fem_gostat_BP_hypo),
    nrow(DMS_24h_fem_gostat_BP_hypo), nrow(DMS_72h_fem_gostat_BP_hypo))

min(nrow(DMS_05h_fem_gostat_BP_hyper), nrow(DMS_1h_fem_gostat_BP_hyper), nrow(DMS_4h_fem_gostat_BP_hyper),
    nrow(DMS_24h_fem_gostat_BP_hyper), nrow(DMS_72h_fem_gostat_BP_hyper), 
    nrow(DMS_05h_fem_gostat_BP_hypo), nrow(DMS_1h_fem_gostat_BP_hypo), nrow(DMS_4h_fem_gostat_BP_hypo),
    nrow(DMS_24h_fem_gostat_BP_hypo), nrow(DMS_72h_fem_gostat_BP_hypo))

max(nrow(DMR_05h_mal_gostat_BP_hyper), nrow(DMR_1h_mal_gostat_BP_hyper), nrow(DMR_4h_mal_gostat_BP_hyper),
    nrow(DMR_24h_mal_gostat_BP_hyper), nrow(DMR_72h_mal_gostat_BP_hyper), 
    nrow(DMR_05h_mal_gostat_BP_hypo), nrow(DMR_1h_mal_gostat_BP_hypo), nrow(DMR_4h_mal_gostat_BP_hypo),
    nrow(DMR_24h_mal_gostat_BP_hypo), nrow(DMR_72h_mal_gostat_BP_hypo))

min(nrow(DMR_05h_mal_gostat_BP_hyper), nrow(DMR_1h_mal_gostat_BP_hyper), nrow(DMR_4h_mal_gostat_BP_hyper),
    nrow(DMR_24h_mal_gostat_BP_hyper), nrow(DMR_72h_mal_gostat_BP_hyper), 
    nrow(DMR_05h_mal_gostat_BP_hypo), nrow(DMR_1h_mal_gostat_BP_hypo), nrow(DMR_4h_mal_gostat_BP_hypo),
    nrow(DMR_24h_mal_gostat_BP_hypo), nrow(DMR_72h_mal_gostat_BP_hypo))

max(nrow(DMS_05h_mal_gostat_BP_hyper), nrow(DMS_1h_mal_gostat_BP_hyper), nrow(DMS_4h_mal_gostat_BP_hyper),
    nrow(DMS_24h_mal_gostat_BP_hyper), nrow(DMS_72h_mal_gostat_BP_hyper), 
    nrow(DMS_05h_mal_gostat_BP_hypo), nrow(DMS_1h_mal_gostat_BP_hypo), nrow(DMS_4h_mal_gostat_BP_hypo),
    nrow(DMS_24h_mal_gostat_BP_hypo), nrow(DMS_72h_mal_gostat_BP_hypo))

min(nrow(DMS_05h_mal_gostat_BP_hyper), nrow(DMS_1h_mal_gostat_BP_hyper), nrow(DMS_4h_mal_gostat_BP_hyper),
    nrow(DMS_24h_mal_gostat_BP_hyper), nrow(DMS_72h_mal_gostat_BP_hyper), 
    nrow(DMS_05h_mal_gostat_BP_hypo), nrow(DMS_1h_mal_gostat_BP_hypo), nrow(DMS_4h_mal_gostat_BP_hypo),
    nrow(DMS_24h_mal_gostat_BP_hypo), nrow(DMS_72h_mal_gostat_BP_hypo))

