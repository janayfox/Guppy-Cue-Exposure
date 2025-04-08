#######################################################
### Goal: Plot and analyze annotation stats 
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")

#load packages
library(ggplot2)
library(dplyr)
library(genomation)
library(DescTools)
library(tibble)
library(ggpubr)

#write functions 
#create function that runs Gtests and post hoc tests
gtest.funct <- function(DMS.data, cpg.data){
  #run initial G test
  print("Initial G test")
  init.test <- GTest(DMS.data, p = cpg.data/(cpg.data[1]+cpg.data[2]+cpg.data[3]+cpg.data[4]))
  print(init.test)
  
  #run post hoc tests
  print("Promoters")
  prom.test <- GTest(c(DMS.data[1], (DMS.data[2]+DMS.data[3]+DMS.data[4])), 
                     p = c(cpg.data[1], (cpg.data[2]+cpg.data[3]+cpg.data[4]))/
                       (cpg.data[1]+cpg.data[2]+cpg.data[3]+cpg.data[4]))
  print(prom.test)
  
  print("Exons")
  ex.test <- GTest(c(DMS.data[2], (DMS.data[1]+DMS.data[3]+DMS.data[4])), 
                   p = c(cpg.data[2], (cpg.data[1]+cpg.data[3]+cpg.data[4]))/
                     (cpg.data[1]+cpg.data[2]+cpg.data[3]+cpg.data[4]))
  print(ex.test)
  
  print("Introns")
  int.test <- GTest(c(DMS.data[3], (DMS.data[1]+DMS.data[2]+DMS.data[4])), 
                    p = c(cpg.data[3], (cpg.data[1]+cpg.data[2]+cpg.data[4]))/
                      (cpg.data[1]+cpg.data[2]+cpg.data[3]+cpg.data[4]))
  print(int.test)
  
  print("Intergenic")
  inter.test <- GTest(c(DMS.data[4], (DMS.data[1]+DMS.data[2]+DMS.data[3])), 
                      p = c(cpg.data[4], (cpg.data[1]+cpg.data[2]+cpg.data[3]))/
                        (cpg.data[1]+cpg.data[2]+cpg.data[3]+cpg.data[4]))
  print(inter.test)
  
  #adjsut p values 
  gtest.p <- c(prom.test$p.value,ex.test$p.value,int.test$p.value,inter.test$p.value)
  p.adjust(gtest.p, method = "hommel")
}

#function for plotting annotation distributions
plot.dist <- function(DMS.data, DMR.data, cpg.data, plotname){
  #format data for plotting 
  DMS.plot.data <- data.frame(DMS.data)
  DMS.plot.data$dist <- "DMS"
  colnames(DMS.plot.data)[1] <- "count"
  DMS.plot.data <- rownames_to_column(DMS.plot.data, "feature")
  
  DMR.plot.data <- data.frame(DMR.data)
  DMR.plot.data$dist <- "DMR"
  colnames(DMR.plot.data)[1] <- "count"
  DMR.plot.data <- rownames_to_column(DMR.plot.data, "feature")
  
  cpg.plot.data <- data.frame(cpg.data)
  cpg.plot.data$dist <- "CpG"
  colnames(cpg.plot.data)[1] <- "count"
  cpg.plot.data <- rownames_to_column(cpg.plot.data, "feature")
  
  plot.data <- rbind(DMS.plot.data, DMR.plot.data, cpg.plot.data)
  
  p <- ggplot(plot.data, aes(x = dist, y = count, fill = feature)) +
    geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) + theme_bw() + xlab(NULL) +
    ylab("Percent") + scale_fill_manual(values = c("exon" = "#B8DE29FF", "intron" = "#1F968BFF", 
                                                   "intergenic" = "#39568CFF", "promoter" = "#481567FF"),
                                        name = "Feature") + 
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) 
  
  ggsave(filename = plotname, plot = p, width = 4, height = 5, units = "in", dpi = 300)
  
  return(p)
}

#format annotation data for plotting 
format.data <- function(DMS.data, DMR.data, cpg.data, time_point){
  DMS.plot.data <- data.frame(DMS.data)
  DMS.plot.data$dist <- "DMS"
  colnames(DMS.plot.data)[1] <- "count"
  DMS.plot.data <- rownames_to_column(DMS.plot.data, "feature")
  
  DMR.plot.data <- data.frame(DMR.data)
  DMR.plot.data$dist <- "DMR"
  colnames(DMR.plot.data)[1] <- "count"
  DMR.plot.data <- rownames_to_column(DMR.plot.data, "feature")
  
  cpg.plot.data <- data.frame(cpg.data)
  cpg.plot.data$dist <- "CpG"
  colnames(cpg.plot.data)[1] <- "count"
  cpg.plot.data <- rownames_to_column(cpg.plot.data, "feature")
  
  bind.data <- rbind(DMS.plot.data, DMR.plot.data, cpg.plot.data)
  bind.data$time_point <- time_point
  return(bind.data)
}

#function for plotting annotation distribution in second way
plot.v2 <- function(DMS.data.05h, DMR.data.05h, cpg.data.05h,
                    DMS.data.1h, DMR.data.1h, cpg.data.1h, 
                    DMS.data.4h, DMR.data.4h, cpg.data.4h,
                    DMS.data.24h, DMR.data.24h, cpg.data.24h,
                    DMS.data.72h, DMR.data.72h, cpg.data.72h){
  #format for plotting 
  data.05h <- format.data(DMS.data.05h, DMR.data.05h, cpg.data.05h, "0.5h")
  data.1h <- format.data(DMS.data.1h, DMR.data.1h, cpg.data.1h, "1h")
  data.4h <- format.data(DMS.data.4h, DMR.data.4h, cpg.data.4h, "4h")
  data.24h <- format.data(DMS.data.24h, DMR.data.24h, cpg.data.24h, "24h")
  data.72h <- format.data(DMS.data.72h, DMR.data.72h, cpg.data.72h, "72h")
  
  plot.data <- rbind(data.05h, data.1h, data.4h, data.24h, data.72h)
  
  #plot 
  p <- ggplot(plot.data, aes(x = dist, y = count, fill = feature)) +
    geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) + theme_bw() + xlab(NULL) +
    ylab("Percent") + scale_fill_manual(values = c("exon" = "#B8DE29FF", "intron" = "#1F968BFF", 
                                                   "intergenic" = "#39568CFF", "promoter" = "#481567FF"),
                                        name = "Feature") + 
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) + facet_grid(~time_point, switch = "x") +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = "white"),
          panel.spacing = unit(-.01,"cm"),
          strip.text.x = element_text(size = 14))
  
  return(p)
}

#load in data
DMS_anno_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_all/DMS_annStats_num_all_fem.RDS")
DMS_anno_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_05h/DMS_annStats_num_05h_fem.RDS")
DMS_anno_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_1h/DMS_annStats_num_1h_fem.RDS")
DMS_anno_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_4h/DMS_annStats_num_4h_fem.RDS")
DMS_anno_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_24h/DMS_annStats_num_24h_fem.RDS")
DMS_anno_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_72h/DMS_annStats_num_72h_fem.RDS")

DMS_anno_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_all/DMS_annStats_num_all_mal.RDS")
DMS_anno_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_05h/DMS_annStats_num_05h_mal.RDS")
DMS_anno_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_1h/DMS_annStats_num_1h_mal.RDS")
DMS_anno_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_4h/DMS_annStats_num_4h_mal.RDS")
DMS_anno_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_24h/DMS_annStats_num_24h_mal.RDS")
DMS_anno_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_72h/DMS_annStats_num_72h_mal.RDS")

DMR_anno_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_all/DMR_annStats_num_all_fem.RDS")
DMR_anno_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_05h/DMR_annStats_num_05h_fem.RDS")
DMR_anno_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_1h/DMR_annStats_num_1h_fem.RDS")
DMR_anno_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_4h/DMR_annStats_num_4h_fem.RDS")
DMR_anno_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_24h/DMR_annStats_num_24h_fem.RDS")
DMR_anno_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_72h/DMR_annStats_num_72h_fem.RDS")

DMR_anno_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_all/DMR_annStats_num_all_mal.RDS")
DMR_anno_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_05h/DMR_annStats_num_05h_mal.RDS")
DMR_anno_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_1h/DMR_annStats_num_1h_mal.RDS")
DMR_anno_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_4h/DMR_annStats_num_4h_mal.RDS")
DMR_anno_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_24h/DMR_annStats_num_24h_mal.RDS")
DMR_anno_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_72h/DMR_annStats_num_72h_mal.RDS")

CpG_anno_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_all/CpG_annStats_num_all_fem.RDS")
CpG_anno_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_05h/CpG_annStats_num_05h_fem.RDS")
CpG_anno_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_1h/CpG_annStats_num_1h_fem.RDS")
CpG_anno_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_4h/CpG_annStats_num_4h_fem.RDS")
CpG_anno_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_24h/CpG_annStats_num_24h_fem.RDS")
CpG_anno_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_72h/CpG_annStats_num_72h_fem.RDS")

CpG_anno_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_all/CpG_annStats_num_all_mal.RDS")
CpG_anno_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_05h/CpG_annStats_num_05h_mal.RDS")
CpG_anno_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_1h/CpG_annStats_num_1h_mal.RDS")
CpG_anno_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_4h/CpG_annStats_num_4h_mal.RDS")
CpG_anno_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_24h/CpG_annStats_num_24h_mal.RDS")
CpG_anno_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res_od/anno_res_72h/CpG_annStats_num_72h_mal.RDS")

## Analyze annotation stats ##
#run g tests
gtest.funct(DMS_anno_all_fem, CpG_anno_all_fem)
gtest.funct(DMS_anno_05h_fem, CpG_anno_05h_fem)
gtest.funct(DMS_anno_1h_fem, CpG_anno_1h_fem)
gtest.funct(DMS_anno_4h_fem, CpG_anno_4h_fem)
gtest.funct(DMS_anno_24h_fem, CpG_anno_24h_fem)
gtest.funct(DMS_anno_72h_fem, CpG_anno_72h_fem)

gtest.funct(DMS_anno_all_mal, CpG_anno_all_mal)
gtest.funct(DMS_anno_05h_mal, CpG_anno_05h_mal)
gtest.funct(DMS_anno_1h_mal, CpG_anno_1h_mal)
gtest.funct(DMS_anno_4h_mal, CpG_anno_4h_mal)
gtest.funct(DMS_anno_24h_mal, CpG_anno_24h_mal)
gtest.funct(DMS_anno_72h_mal, CpG_anno_72h_mal)

gtest.funct(DMR_anno_all_fem, CpG_anno_all_fem)
gtest.funct(DMR_anno_05h_fem, CpG_anno_05h_fem)
gtest.funct(DMR_anno_1h_fem, CpG_anno_1h_fem)
gtest.funct(DMR_anno_4h_fem, CpG_anno_4h_fem)
gtest.funct(DMR_anno_24h_fem, CpG_anno_24h_fem)
gtest.funct(DMR_anno_72h_fem, CpG_anno_72h_fem)

gtest.funct(DMR_anno_all_mal, CpG_anno_all_mal)
gtest.funct(DMR_anno_05h_mal, CpG_anno_05h_mal)
gtest.funct(DMR_anno_1h_mal, CpG_anno_1h_mal)
gtest.funct(DMR_anno_4h_mal, CpG_anno_4h_mal)
gtest.funct(DMR_anno_24h_mal, CpG_anno_24h_mal)
gtest.funct(DMR_anno_72h_mal, CpG_anno_72h_mal)

#plot annotation distribution 
plot.dist(DMS_anno_all_fem, DMR_anno_all_fem, CpG_anno_all_fem,"./shortTerm_exp/plots/finalized_tiff/anno_barplots_od/DMS_allTP_fem_anno_barplot.tiff")
fem_05h_plot <- plot.dist(DMS_anno_05h_fem, DMR_anno_05h_fem, CpG_anno_05h_fem,"./shortTerm_exp/plots/finalized_tiff/anno_barplots/DMS_05h_fem_anno_barplot.tiff")
fem_1h_plot <- plot.dist(DMS_anno_1h_fem, DMR_anno_1h_fem, CpG_anno_1h_fem,"./shortTerm_exp/plots/finalized_tiff/anno_barplots/DMS_1h_fem_anno_barplot.tiff")
fem_4h_plot <- plot.dist(DMS_anno_4h_fem, DMR_anno_4h_fem, CpG_anno_4h_fem,"./shortTerm_exp/plots/finalized_tiff/anno_barplots/DMS_4h_fem_anno_barplot.tiff")
fem_24h_plot <- plot.dist(DMS_anno_24h_fem, DMR_anno_24h_fem, CpG_anno_24h_fem,"./shortTerm_exp/plots/finalized_tiff/anno_barplots/DMS_24h_fem_anno_barplot.tiff")
fem_72h_plot <- plot.dist(DMS_anno_72h_fem, DMR_anno_72h_fem, CpG_anno_72h_fem,"./shortTerm_exp/plots/finalized_tiff/anno_barplots/DMS_72h_fem_anno_barplot.tiff")

plot.dist(DMS_anno_all_mal, DMR_anno_all_mal, CpG_anno_all_mal,"./shortTerm_exp/plots/finalized_tiff/anno_barplots_od/DMS_allTP_mal_anno_barplot.tiff")
mal_05h_plot <- plot.dist(DMS_anno_05h_mal, DMR_anno_05h_mal, CpG_anno_05h_mal,"./shortTerm_exp/plots/finalized_tiff/anno_barplots/DMS_05h_mal_anno_barplot.tiff")
mal_1h_plot <- plot.dist(DMS_anno_1h_mal, DMR_anno_1h_mal, CpG_anno_1h_mal,"./shortTerm_exp/plots/finalized_tiff/anno_barplots/DMS_1h_mal_anno_barplot.tiff")
mal_4h_plot <- plot.dist(DMS_anno_4h_mal, DMR_anno_4h_mal, CpG_anno_4h_mal,"./shortTerm_exp/plots/finalized_tiff/anno_barplots/DMS_4h_mal_anno_barplot.tiff")
mal_24h_plot <- plot.dist(DMS_anno_24h_mal, DMR_anno_24h_mal, CpG_anno_24h_mal,"./shortTerm_exp/plots/finalized_tiff/anno_barplots/DMS_24h_mal_anno_barplot.tiff")
mal_72h_plot <- plot.dist(DMS_anno_72h_mal, DMR_anno_72h_mal, CpG_anno_72h_mal,"./shortTerm_exp/plots/finalized_tiff/anno_barplots/DMS_72h_mal_anno_barplot.tiff")

#make panels
panel <- ggarrange(fem_05h_plot, fem_1h_plot, fem_4h_plot, fem_24h_plot, fem_72h_plot,
                   mal_05h_plot, mal_1h_plot, mal_4h_plot, mal_24h_plot, mal_72h_plot,
                           common.legend = TRUE, nrow = 2, ncol = 5, legend = "bottom",
                           labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"))

ggsave(filename = "./shortTerm_exp/plots/finalized_tiff/anno_barplots/all_panel.tiff", plot = panel, width = 12, height = 6, units = "in", dpi = 300)

#try plotting in a different way 
fem.plot <- plot.v2(DMS_anno_05h_fem, DMR_anno_05h_fem, CpG_anno_05h_fem,
                    DMS_anno_1h_fem, DMR_anno_1h_fem, CpG_anno_1h_fem,
                    DMS_anno_4h_fem, DMR_anno_4h_fem, CpG_anno_4h_fem,
                    DMS_anno_24h_fem, DMR_anno_24h_fem, CpG_anno_24h_fem,
                    DMS_anno_72h_fem, DMR_anno_72h_fem, CpG_anno_72h_fem)

mal.plot <- plot.v2(DMS_anno_05h_mal, DMR_anno_05h_mal, CpG_anno_05h_mal,
                    DMS_anno_1h_mal, DMR_anno_1h_mal, CpG_anno_1h_mal,
                    DMS_anno_4h_mal, DMR_anno_4h_mal, CpG_anno_4h_mal,
                    DMS_anno_24h_mal, DMR_anno_24h_mal, CpG_anno_24h_mal,
                    DMS_anno_72h_mal, DMR_anno_72h_mal, CpG_anno_72h_mal)

panel <- ggarrange(fem.plot, mal.plot,
                         labels = c("A", "B"), nrow = 2, ncol = 1,
                   common.legend = TRUE, legend = "bottom")


ggsave(filename = "./shortTerm_exp/plots/finalized_tiff/anno_barplots/panel_anno_plot.tiff", plot = panel, width = 8, height = 6, units = "in", dpi = 300)
