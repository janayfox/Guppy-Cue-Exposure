#######################################################
### Goal: Plot heatmaps and do cluster analysis
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("ggpubr")
#install.packages("ggupset")
#install.packages("VennDiagram")

#load packaes
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggupset)
library(VennDiagram)
library(methylKit)
library(pheatmap)

#create function that makes a column with a name for each DMS/DMR
nameSite <- function(data){
  data$site_name <- paste(data$chr, data$start, sep = "_")
  return(data)
}

#function that generates a percent methylation matrix
get_percMeth_matrix <- function(diffMethfile, methFile){
  
  #read data files
  diffMeth <- readRDS(diffMethfile)
  meth <- readRDS(methFile)
  
  #calculate percent methylation 
  perc_meth <- as.data.frame(percMethylation(meth))
  
  meth <- nameSite(meth)
  diffMeth <- nameSite(diffMeth)
  
  #add rownames to perc meth matrix  
  rownames(perc_meth) <- meth$site_name
  
  #subset only DMS from meth matrix 
  perc_meth <- perc_meth[rownames(perc_meth) %in% diffMeth$site_name, ]
  
  return(perc_meth)
  
}

#run for all datasets
DMS_fem_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSdiffmeth_all_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSmeth_all_fem_5X.RDS")

DMS_fem_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_05h/DMSdiffmeth_05h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_05h/DMSmeth_05h_fem_5X.RDS")

DMS_fem_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_1h/DMSdiffmeth_1h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_1h/DMSmeth_1h_fem_5X.RDS")

DMS_fem_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_4h/DMSdiffmeth_4h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_4h/DMSmeth_4h_fem_5X.RDS")

DMS_fem_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_24h/DMSdiffmeth_24h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_24h/DMSmeth_24h_fem_5X.RDS")

DMS_fem_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_72h/DMSdiffmeth_72h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_72h/DMSmeth_72h_fem_5X.RDS")

DMS_mal_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSdiffmeth_all_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSmeth_all_mal_5X.RDS")

DMS_mal_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_05h/DMSdiffmeth_05h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMS_res_05h/DMSmeth_05h_mal_5X.RDS")

DMS_mal_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_1h/DMSdiffmeth_1h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_1h/DMSmeth_1h_mal_5X.RDS")

DMS_mal_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_4h/DMSdiffmeth_4h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMS_res_4h/DMSmeth_4h_mal_5X.RDS")

DMS_mal_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_24h/DMSdiffmeth_24h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMS_res_24h/DMSmeth_24h_mal_5X.RDS")

DMS_mal_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMS_res_72h/DMSdiffmeth_72h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMS_res_72h/DMSmeth_72h_mal_5X.RDS")

DMR_fem_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_all/DMRdiffmeth_all_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_all/DMRmeth_all_fem_5X.RDS")

DMR_fem_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_05h/DMRdiffmeth_05h_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_05h/DMRmeth_05h_fem_5X.RDS")

DMR_fem_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_1h/DMRdiffmeth_1h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMR_res_1h/DMRmeth_1h_fem_5X.RDS")

DMR_fem_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_4h/DMRdiffmeth_4h_fem_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMR_res_4h/DMRmeth_4h_fem_5X.RDS")

DMR_fem_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_24h/DMRdiffmeth_24h_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_24h/DMRmeth_24h_fem_5X.RDS")

DMR_fem_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_72h/DMRdiffmeth_72h_fem_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_72h/DMRmeth_72h_fem_5X.RDS")

DMR_mal_all_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_all/DMRdiffmeth_all_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_all/DMRmeth_all_mal_5X.RDS")

DMR_mal_05h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_05h/DMRdiffmeth_05h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_05h/DMRmeth_05h_mal_5X.RDS")

DMR_mal_1h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_1h/DMRdiffmeth_1h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMR_res_1h/DMRmeth_1h_mal_5X.RDS")

DMR_mal_4h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_4h/DMRdiffmeth_4h_mal_5X.RDS",
                                           "./shortTerm_exp/data/methylkit_res/DMR_res_4h/DMRmeth_4h_mal_5X.RDS")

DMR_mal_24h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_24h/DMRdiffmeth_24h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_24h/DMRmeth_24h_mal_5X.RDS")

DMR_mal_72h_percMeth <- get_percMeth_matrix("./shortTerm_exp/data/methylkit_res/DMR_res_72h/DMRdiffmeth_72h_mal_5X.RDS",
                                            "./shortTerm_exp/data/methylkit_res/DMR_res_72h/DMRmeth_72h_mal_5X.RDS")

#get methylation matrix for all sites 
#read data files
meth_allsites <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSmeth_all_5X.RDS")

#calculate percent methylation 
perc_meth_allsites <- as.data.frame(percMethylation(meth_allsites))

#add site names
meth_allsites <- nameSite(meth_allsites)

#add rownames to perc meth matrix  
rownames(perc_meth_allsites) <- meth_allsites$site_name

#remove extra file 
rm(meth_allsites)

##Run PCAs## 
#get lists of all DMSs for each sex
DMS_combined <- c(row.names(DMS_fem_05h_percMeth), row.names(DMS_fem_1h_percMeth), row.names(DMS_fem_4h_percMeth), 
                  row.names(DMS_fem_24h_percMeth), row.names(DMS_fem_72h_percMeth),
                  row.names(DMS_mal_05h_percMeth), row.names(DMS_mal_1h_percMeth), row.names(DMS_mal_4h_percMeth), 
                  row.names(DMS_mal_24h_percMeth), row.names(DMS_mal_72h_percMeth))

DMS_combined <- unique(DMS_combined) #remove duplicates

DMR_combined <- c(row.names(DMR_fem_05h_percMeth), row.names(DMR_fem_1h_percMeth), row.names(DMR_fem_4h_percMeth), 
                  row.names(DMR_fem_24h_percMeth), row.names(DMR_fem_72h_percMeth),
                  row.names(DMR_mal_05h_percMeth), row.names(DMR_mal_1h_percMeth), row.names(DMR_mal_4h_percMeth), 
                  row.names(DMR_mal_24h_percMeth), row.names(DMR_mal_72h_percMeth))

DMR_combined <- unique(DMR_combined) #remove duplicates

#extract DMSs and DMRs
DMS_combined_meth <- perc_meth_allsites %>% filter(row.names(perc_meth_allsites) %in% DMS_combined)
DMR_combined_meth <- perc_meth_allsites %>% filter(row.names(perc_meth_allsites) %in% DMR_combined)

#make function to transpose data 
prep_data <- function(df){
  #transpose data
  t.data <- t(df)
  rownames(t.data) <- colnames(df)
  colnames(t.data) <- rownames(df)
  return(as.data.frame(t.data))
}

t_DMS_combined_meth <- prep_data(DMS_combined_meth)
t_DMR_combined_meth <- prep_data(DMR_combined_meth)
t_DMS_fem_all <- prep_data(DMS_fem_all_percMeth)
t_DMS_mal_all <- prep_data(DMS_mal_all_percMeth)

#run PCA
pca_DMS_combined <- prcomp(t_DMS_combined_meth, center = TRUE)
pca_DMR_combined <- prcomp(t_DMR_combined_meth, center = TRUE)
pca_DMS_fem_all <- prcomp(t_DMS_fem_all, center = TRUE)
pca_DMS_mal_all <- prcomp(t_DMS_mal_all, center = TRUE)

#view results 
summary(pca_DMS_combined)
summary(pca_DMR_combined)
summary(pca_DMS_fem_all)
summary(pca_DMS_mal_all)

#make lists of IDs for each time point 
IDs_05h <- c(colnames(DMS_fem_05h_percMeth), colnames(DMS_mal_05h_percMeth))
IDs_1h <- c(colnames(DMS_fem_1h_percMeth), colnames(DMS_mal_1h_percMeth))
IDs_4h <- c(colnames(DMS_fem_4h_percMeth), colnames(DMS_mal_4h_percMeth))
IDs_24h <- c(colnames(DMS_fem_24h_percMeth), colnames(DMS_mal_24h_percMeth))
IDs_72h <- c(colnames(DMS_fem_72h_percMeth), colnames(DMS_mal_72h_percMeth))

#add on site and collection data 
add_var <- function(data, pca) {
  #add first four PCs to df
  data <- cbind(data, pca$x[,1:4])
  
  #convert rownames to ID column 
  data.new <- tibble::rownames_to_column(data, "ID")
  
  #add sex
  data.new$sex <- data.new$ID
  data.new$sex <- str_sub(data.new$sex, -1)
  
  #add cue
  data.new$cue <- data.new$ID
  data.new$cue <- str_sub(data.new$cue, 4, 4)
  
  #add tank
  data.new$tank <- data.new$ID
  data.new$tank <- str_sub(data.new$tank, 4, 7)
  data.new$tank <- gsub("F", "", data.new$tank)
  data.new$tank <- gsub("M", "", data.new$tank)
  
  #add time point 
  data.new$time_point <- "NA"
  data.new <- within(data.new, time_point[ID %in% IDs_05h] <- "05h")
  data.new <- within(data.new, time_point[ID %in% IDs_1h] <- "1h")
  data.new <- within(data.new, time_point[ID %in% IDs_4h] <- "4h")
  data.new <- within(data.new, time_point[ID %in% IDs_24h] <- "24h")
  data.new <- within(data.new, time_point[ID %in% IDs_72h] <- "72h")
  
  #put rownames back on 
  rownames(data.new) <- rownames(data)
  return(data.new)
}

t_DMS_combined_meth <- add_var(t_DMS_combined_meth, pca_DMS_combined)
t_DMR_combined_meth <- add_var(t_DMR_combined_meth, pca_DMR_combined)
t_DMS_fem_all <- add_var(t_DMS_fem_all, pca_DMS_fem_all)
t_DMS_mal_all <- add_var(t_DMS_mal_all, pca_DMS_mal_all)

#plot PCAs
plot_pca <- function(pca.data, pca.x, pca.y, xlabel, ylabel){
  ggplot(data = pca.data, aes(x = pca.x, y = pca.y, color = cue, shape = time_point)) + theme_bw() + geom_point() +
    scale_color_manual(name = "Cue",
                       values = c("#E14D2A", "#31C6D4")) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    xlab(xlabel) + ylab(ylabel) +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16), legend.text = element_text(size=18), 
          legend.title = element_text(size = 18))
}

plot_pca(t_DMS_combined_meth, t_DMS_combined_meth$PC1, t_DMS_combined_meth$PC2, "PC1 (26.31%)", "PC2 (13.46%)")

##Make heat maps ## 
#cluster samples
sample.clust <- hclust(dist(t(perc_meth), method = "euclidean"), method = "ward.D2")

#generate color pallette 
cor.color <- colorRampPalette(c("#FFFF00", "black", "#0000FF"))(60) 

#calculate breaks 
breaks = c(seq(-6, 0, length.out = 30),seq(0.01, 6, length.out = 30))

set.breaks <- function(expr.data){
  new.breaks <- breaks
  new.breaks[length(breaks)] <- max(max(expr.data),max(breaks))
  new.breaks[1] <- min(min(expr.data),min(breaks))
  return(new.breaks)
}


#make dataframe for sample annotations 
sample.info <- data.frame(sample = colnames(perc_meth))
sample.info[substr(sample.info$sample,4,4)=="A", "Cue"] <- "Alarm Cue"  
sample.info[substr(sample.info$sample,4,4)=="C", "Cue"] <- "Control Cue"  

rownames(sample.info) <- sample.info[,1]
sample.info[,1] <- NULL

anno.colors <- list(Cue = c("Alarm Cue" = "#E14D2A", "Control Cue" = "#31C6D4"))

#maybe would be good to do gene clustering again in a similar way?? can look at them over time 

# to put genes in order of clustering without having dendrogram, 
#I can do separate clustering of the genes however I want, add the cluster information on,
#then arrange the rows by the cluster  info, then it will plot them in order that the rows are in

#plot heatmaps

p <-  pheatmap(perc_meth, cluster_rows = TRUE, cluster_cols = sample.clust, show_colnames = FALSE, show_rownames = FALSE, scale = "row",
               color = cor.color, annotation_col = sample.info, annotation_colors = anno.colors, treeheight_row = 0)
p
           
color = cor.color, breaks = breaks.samp, annotation_row = gene.anno,
           annotation_col = sample.anno, annotation_colors = anno.colors,
           annotation_names_col = FALSE, annotation_names_row = FALSE, legend_breaks = c(-6,-3,0,3,6), border_color = NA)

#set annotation colours

