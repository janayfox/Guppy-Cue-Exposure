#######################################################
### Goal: Plot venn diagram and heatmaps and do cluster analysis
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("ggpubr")
#install.packages("VennDiagram")

#load packaes
library(ggplot2)
library(dplyr)
library(ggpubr)
library(VennDiagram)
library(methylKit)
library(pheatmap)
library(stringr)

#read data 
DMS_diffmeth_fem <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSdiffmeth_fem_5X_data.RDS")
DMS_diffmeth_mal <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSdiffmeth_mal_5X_data.RDS")
DMS_diffmeth_all <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSdiffmeth_all_5X_data.RDS")

DMR_diffmeth_fem <- readRDS("./dev_exp/data/methylkit_res/DMR_res/DMRdiffmeth_fem_5X_data.RDS")
DMR_diffmeth_mal <- readRDS("./dev_exp/data/methylkit_res/DMR_res/DMRdiffmeth_mal_5X_data.RDS")
DMR_diffmeth_all <- readRDS("./dev_exp/data/methylkit_res/DMR_res/DMRdiffmeth_all_5X_data.RDS")

## Setup functions ##
# function that makes a column with a name for each DMS/DMR
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

#function that plots hyper vs hypo plots 
plot_hyp <- function(data, ylab){
  p <- ggplot(data, aes(x = type, y = number, fill = direction)) + theme_bw() +
    geom_bar(stat = "identity", position = "dodge", color = "black") + 
    labs(x = "Comparison", y = ylab, fill = "Direction of Methylation") + scale_fill_manual(values = c("lightgreen", "orange")) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text = element_text(size = 11),
          legend.title = element_text(size = 12)) 
}

plot_hyp_all <- function(data){
  p <- ggplot(data, aes(x = type, y = number, fill = Direction)) + theme_bw() +
    geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Percent") + scale_fill_manual(values = c("lightgreen", "orange")) + 
    theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16), legend.text = element_text(size = 14),
          legend.title = element_text(size = 15))   
}

#function that plots manhattan plots
plotManhat <- function(df, plotname){
  #set colors
  color_val <- rep(c("lightgreen", "skyblue"), length.out = 23)
  
  df$chr <- as.character(df$chr)
  
  #convert chromosome names for plotting
  df$chr[df$chr == "NC_024331.1"] <- 1
  df$chr[df$chr == "NC_024332.1"] <- 2
  df$chr[df$chr == "NC_024333.1"] <- 3
  df$chr[df$chr == "NC_024334.1"] <- 4
  df$chr[df$chr == "NC_024335.1"] <- 5
  df$chr[df$chr == "NC_024336.1"] <- 6
  df$chr[df$chr == "NC_024337.1"] <- 7
  df$chr[df$chr == "NC_024338.1"] <- 8
  df$chr[df$chr == "NC_024339.1"] <- 9
  df$chr[df$chr == "NC_024340.1"] <- 10
  df$chr[df$chr == "NC_024341.1"] <- 11
  df$chr[df$chr == "NC_024342.1"] <- 12
  df$chr[df$chr == "NC_024343.1"] <- 13
  df$chr[df$chr == "NC_024344.1"] <- 14
  df$chr[df$chr == "NC_024345.1"] <- 15
  df$chr[df$chr == "NC_024346.1"] <- 16
  df$chr[df$chr == "NC_024347.1"] <- 17
  df$chr[df$chr == "NC_024348.1"] <- 18
  df$chr[df$chr == "NC_024349.1"] <- 19
  df$chr[df$chr == "NC_024350.1"] <- 20
  df$chr[df$chr == "NC_024351.1"] <- 21
  df$chr[df$chr == "NC_024352.1"] <- 22
  df$chr[df$chr == "NC_024353.1"] <- 23
  
  #convert chr column to numeric
  df$chr <- as.numeric(df$chr)
  
  #calculate cumulative position of SNP
  df.don <- df %>%
    #compute chr size
    group_by(chr) %>%
    summarise(chr_len=max(start)) %>%
    #calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    #add this info to the initial dataset
    left_join(df, ., by=c("chr"="chr")) %>%
    #add a cumulative position of each SNP
    arrange(chr, start) %>%
    mutate(DMScum=start+tot)
  
  print(df.don)
  
  #make axis
  axis <- df.don %>% group_by(chr) %>% summarize(center=(max(DMScum) + min(DMScum))/2)
  
  #make Manhattan plots
  p <- ggplot(df.don, aes(x=DMScum, y=meth.diff)) +
    geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = color_val) +
    scale_x_continuous(label = axis$chr, breaks = axis$center) +
    theme_light() + ylab("Change in % methylation") + xlab("Chromosome") +
    geom_hline(yintercept = 0, color = "red", linetype="dashed") +
    theme(
      legend.position="none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  ggsave(filename = plotname, plot = p, width = 7, height = 4, units = "in", dpi = 300)
  
  return(p)
}

#function for adding on tank, cue and sex data 
add_var <- function(data, pca) {
  #make sure data is dataframe 
  data <- as.data.frame(data)
  
  #add first four PCs to df
  data <- cbind(data, pca$x[,1:4])
  
  #convert rownames to ID column 
  data.new <- tibble::rownames_to_column(data, "ID")
  
  #add sex
  data.new$sex <- data.new$ID
  data.new$sex <- str_sub(data.new$sex, -2, -2)
  
  #add cue
  data.new$cue <- data.new$ID
  data.new$cue <- str_sub(data.new$cue, 2, 2)
  
  #add tank
  data.new$tank <- data.new$ID
  data.new$tank <- str_sub(data.new$tank, 2, 4)
  data.new$tank <- gsub("F", "", data.new$tank)
  data.new$tank <- gsub("M", "", data.new$tank)
  
  #put rownames back on 
  rownames(data.new) <- rownames(data)
  return(data.new)
}


#make function to transpose data 
prep_data <- function(df){
  #remove NAs
  df <- na.omit(df)
  
  #transpose data
  t.df <- t(df)
  
  return(t.df)
}

#function for plotting PCAs 
plot_pca <- function(pca.data, pca.x, pca.y, xlabel, ylabel, plotname){
  p <- ggplot(data = pca.data, aes(x = pca.x, y = pca.y, color = cue, shape = tank, group = cue)) + theme_bw() + geom_point() +
    scale_color_manual(name = "Cue",
                       values = c("#E14D2A", "lightblue"), labels = c("A" = "Alarm Cue", "C" = "Control")) +
    scale_shape_manual(values = c(15, 25, 17, 18, 19, 8, 15, 25, 17, 18, 19, 8)) +
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    xlab(xlabel) + ylab(ylabel) + labs(shape = "Tank") + stat_ellipse() +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16), legend.text = element_text(size=14),
          legend.title = element_text(size = 18))
  
  ggsave(filename = plotname, plot = p, width = 6, height = 5, units = "in", dpi = 300)
  
  return(p)
}

plot_pca_V2 <- function(pca.data, pca.x, pca.y, xlabel, ylabel, plotname){
  p <- ggplot(data = pca.data, aes(x = pca.x, y = pca.y, color = cue, shape = sex, group = cue)) + theme_bw() + geom_point() +
    scale_color_manual(name = "Cue",
                       values = c("#E14D2A", "lightblue"), labels = c("A" = "Alarm Cue", "C" = "Control")) +
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    xlab(xlabel) + ylab(ylabel) + labs(shape = "Sex") + stat_ellipse() +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16), legend.text = element_text(size=14),
          legend.title = element_text(size = 18))
  
  ggsave(filename = plotname, plot = p, width = 6, height = 5, units = "in", dpi = 300)
}

#make function to plot heatmaps of individual time points
plot_heatmap <- function(data, plotname){
  
  #replace NAs with 0
  data[is.na(data)] <- 0
  
  #cluster samples
  sample.clust <- hclust(dist(t(data), method = "euclidean"), method = "ward.D2")
  
  #generate color pallette for heatmap
  cor.color <- colorRampPalette(c("#FFFF00", "black", "#0000FF"))(100)
  
  #make dataframe for sample annotations
  sample.info <- data.frame(sample = colnames(data))
  sample.info[str_sub(sample.info$sample, 2, 2)=="A", "Cue"] <- "Alarm Cue"
  sample.info[str_sub(sample.info$sample, 2 ,2)=="C", "Cue"] <- "Control"
  # sample.info$Tank <- "NA"
  # sample.info <- within(sample.info, Tank[sample %in% AC2_list] <- "AC2")
  # sample.info <- within(sample.info, Tank[sample %in% AC3_list] <- "AC3")
  # sample.info <- within(sample.info, Tank[sample %in% AC4_list] <- "AC4")
  # sample.info <- within(sample.info, Tank[sample %in% AC5_list] <- "AC5")
  # sample.info <- within(sample.info, Tank[sample %in% AC6_list] <- "AC6")
  # sample.info <- within(sample.info, Tank[sample %in% AC7_list] <- "AC7")
  # sample.info <- within(sample.info, Tank[sample %in% C2_list] <- "C2")
  # sample.info <- within(sample.info, Tank[sample %in% C3_list] <- "C3")
  # sample.info <- within(sample.info, Tank[sample %in% C4_list] <- "C4")
  # sample.info <- within(sample.info, Tank[sample %in% C5_list] <- "C5")
  # sample.info <- within(sample.info, Tank[sample %in% C6_list] <- "C6")
  # sample.info <- within(sample.info, Tank[sample %in% C7_list] <- "C7")
  rownames(sample.info) <- sample.info[,1]
  sample.info[,1] <- NULL
  
  #set annotation colors
  anno.colors <- list(Cue = c("Alarm Cue" = "#E14D2A", "Control" = "lightblue"))
  
                      # Tank = c("C2" = "#03045E", "C3" = "#023E8A", "C4" = "#0077B6",
                      #          "C5" = "#183EFA", "C6" = "#00FFEF", "C7" = "#82EEFD",
                      #          "AC2" = "#541E1B", "AC3" = "#960018", "AC4" = "#E3242B",
                      #          "AC5" = "#FF6347", "AC6" = "#FC7676", "AC7" = "#FFCBCB"))
  
  #scale by rows
  scal.data <- apply(data, 1, scale)
  scal.data <- t(scal.data)
  colnames(scal.data) <- colnames(data)
  
  #calculate breaks 
  breaks = c(seq(-3, 0, length.out = 50),seq(0.001, 3, length.out = 50))
  new.breaks <- breaks
  new.breaks[length(breaks)] <- max(max(scal.data),max(breaks))
  new.breaks[1] <- min(min(scal.data),min(breaks))
  
  #plot heatmap
  p <-  pheatmap(scal.data, cluster_rows = TRUE, cluster_cols = sample.clust,
                 show_colnames = FALSE, show_rownames = FALSE, scale = "none",
                 border_color = NA,
                 breaks = new.breaks, legend_breaks = c(-3,-2,-1,0,1,2,3),
                 color = cor.color, annotation_col = sample.info,
                 annotation_colors = anno.colors, treeheight_row = 0)
  
  #save plot
  tiff(plotname, units="in", width = 6, height = 6, res = 600)
  print(p)
  dev.off()
}

#function for plotting heatmaps of both sexes together
plot_heatmap_all <- function(data, plotname){
  
  #replace NAs with 0
  data[is.na(data)] <- 0
  
  #cluster samples
  sample.clust <- hclust(dist(t(data), method = "euclidean"), method = "ward.D2")
  
  #generate color pallette for heatmap
  cor.color <- colorRampPalette(c("#FFFF00", "black", "#0000FF"))(100)
  
  #make dataframe for sample annotations
  sample.info <- data.frame(sample = colnames(data))
  sample.info[str_sub(sample.info$sample, 2, 2)=="A", "Cue"] <- "Alarm Cue"
  sample.info[str_sub(sample.info$sample, 2 , 2)=="C", "Cue"] <- "Control"
  sample.info[str_sub(sample.info$sample, -2, -2)=="F", "Sex"] <- "Female"
  sample.info[str_sub(sample.info$sample, -2 , -2)=="M", "Sex"] <- "Male"
  rownames(sample.info) <- sample.info[,1]
  sample.info[,1] <- NULL
  
  #set annotation colors
  anno.colors <- list(Cue = c("Alarm Cue" = "#E14D2A", "Control" = "#31C6D4"),
                      Sex = c("Female" = "#FFE17B", "Male" = "skyblue"))
  
  #scale by rows
  scal.data <- apply(data, 1, scale)
  scal.data <- t(scal.data)
  colnames(scal.data) <- colnames(data)
  
  #calculate breaks 
  breaks = c(seq(-3, 0, length.out = 50),seq(0.001, 3, length.out = 50))
  new.breaks <- breaks
  new.breaks[length(breaks)] <- max(max(scal.data),max(breaks))
  new.breaks[1] <- min(min(scal.data),min(breaks))
  
  #plot heatmap
  p <-  pheatmap(scal.data, cluster_rows = TRUE, cluster_cols = sample.clust,
                 show_colnames = FALSE, show_rownames = FALSE, scale = "none",
                 border_color = NA,
                 breaks = new.breaks, legend_breaks = c(-3,-2,-1,0,1,2,3),
                 color = cor.color, annotation_col = sample.info,
                 annotation_colors = anno.colors, treeheight_row = 0)
  
  #save plot
  tiff(plotname, units="in", width = 6, height = 6, res = 600)
  print(p)
  dev.off()
}

## Plot venn diagrams ##
#add site name to diffmeth
DMS_diffmeth_fem <- nameSite(DMS_diffmeth_fem)
DMS_diffmeth_mal <- nameSite(DMS_diffmeth_mal)
DMS_diffmeth_all <- nameSite(DMS_diffmeth_all)

DMR_diffmeth_fem <- nameSite(DMR_diffmeth_fem)
DMR_diffmeth_mal <- nameSite(DMR_diffmeth_mal)
DMR_diffmeth_all <- nameSite(DMR_diffmeth_all)

#make lists
DMS.fem <- DMS_diffmeth_fem$site_name
DMS.mal <- DMS_diffmeth_mal$site_name
DMS.all <- DMS_diffmeth_all$site_name

DMR.fem <- DMR_diffmeth_fem$site_name
DMR.mal <- DMR_diffmeth_mal$site_name
DMR.all <- DMR_diffmeth_all$site_name

DMS.femVmal <- list(DMS.mal, DMS.fem)
DMS.femVmalVall <- list(DMS.mal, DMS.fem, DMS.all)

DMR.femVmal <- list(DMR.mal, DMR.fem)
DMR.femVmalVall <- list(DMR.mal, DMR.fem, DMR.all)
# 
# venn.diagram(DMS.femVmal, category.names = c("Males", "Females"),
#             filename = "./dev_exp/plots/finalized_tiff/venn_diagrams/DMS_femVmal_venn.png",
#             imagetype = "tiff", height = 2000, width = 2000, resolution = 600,
#             lwd = 2, col = c("skyblue", "#FFE17B"),
#             fill = c(alpha("skyblue", 0.8), alpha("#FFE17B", 0.8)),
#             cex = 1, cat.cex = 1, scaled = FALSE, fontfamily = "sans",
#             cat.fontfamily = "sans",cat.default.pos = "outer",cat.pos = c(-27, 27),
#             cat.dist = c(0.055, 0.055))
# 
# venn.diagram(DMS.femVmalVall, category.names = c("Males", "Females", "All"),
#              filename = "./dev_exp/plots/finalized_tiff/venn_diagrams/DMS_femVmalVall_venn.png",
#              imagetype = "tiff", height = 2000, width = 2000, resolution = 600,
#              lwd = 2, col = c("skyblue", "#FFE17B", "lightgreen"),
#              fill = c(alpha("skyblue", 0.8), alpha("#FFE17B", 0.8), alpha("lightgreen", 0.8)),
#              cex = 1, cat.cex = 1, scaled = FALSE, fontfamily = "sans",
#              cat.fontfamily = "sans",cat.default.pos = "outer",
#              cat.dist = c(0.055, 0.055, 0.055))
# 
# venn.diagram(DMR.femVmal, category.names = c("Males", "Females"),
#              filename = "./dev_exp/plots/finalized_tiff/venn_diagrams/DMR_femVmal_venn.png",
#              imagetype = "tiff", height = 2000, width = 2000, resolution = 600,
#              lwd = 2, col = c("skyblue", "#FFE17B"),
#              fill = c(alpha("skyblue", 0.8), alpha("#FFE17B", 0.8)),
#              cex = 1, cat.cex = 1, scaled = FALSE, fontfamily = "sans",
#              cat.fontfamily = "sans",cat.default.pos = "outer",cat.pos = c(-27, 27),
#              cat.dist = c(0.055, 0.055))
# 
# venn.diagram(DMR.femVmalVall, category.names = c("Males", "Females", "All"),
#              filename = "./dev_exp/plots/finalized_tiff/venn_diagrams/DMR_femVmalVall_venn.png",
#              imagetype = "tiff", height = 2000, width = 2000, resolution = 600,
#              lwd = 2, col = c("skyblue", "#FFE17B", "lightgreen"),
#              fill = c(alpha("skyblue", 0.8), alpha("#FFE17B", 0.8), alpha("lightgreen", 0.8)),
#              cex = 1, cat.cex = 1, scaled = FALSE, fontfamily = "sans",
#              cat.fontfamily = "sans",cat.default.pos = "outer", cat.pos = c(-27,27,-30),
#              cat.dist = c(0.055, 0.055, 0.055))

# ## Plot barplots of hyper vs hypo ## 
#check numbers of hyper vs hypo
nrow(subset(DMS_diffmeth_fem, meth.diff > 0))
nrow(subset(DMS_diffmeth_fem, meth.diff < 0))

nrow(subset(DMS_diffmeth_mal, meth.diff > 0))
nrow(subset(DMS_diffmeth_mal, meth.diff < 0))

nrow(subset(DMS_diffmeth_all, meth.diff > 0))
nrow(subset(DMS_diffmeth_all, meth.diff < 0))

nrow(subset(DMR_diffmeth_fem, meth.diff > 0))
nrow(subset(DMR_diffmeth_fem, meth.diff < 0))

nrow(subset(DMR_diffmeth_mal, meth.diff > 0))
nrow(subset(DMR_diffmeth_mal, meth.diff < 0))

nrow(subset(DMR_diffmeth_all, meth.diff > 0))
nrow(subset(DMR_diffmeth_all, meth.diff < 0))

#check max and min methylation changes in DMRs
min(abs(DMS_diffmeth_fem$meth.diff))
max(abs(DMS_diffmeth_fem$meth.diff))
min(abs(DMR_diffmeth_fem$meth.diff))
max(abs(DMR_diffmeth_fem$meth.diff))

min(abs(DMS_diffmeth_mal$meth.diff))
max(abs(DMS_diffmeth_mal$meth.diff))
min(abs(DMR_diffmeth_mal$meth.diff))
max(abs(DMR_diffmeth_mal$meth.diff))

min(abs(DMS_diffmeth_all$meth.diff))
max(abs(DMS_diffmeth_all$meth.diff))
min(abs(DMR_diffmeth_all$meth.diff))
max(abs(DMR_diffmeth_all$meth.diff))

#get dataframes
DMS.hyp.data <- data.frame(direction = c("Hyper", "Hypo", "Hyper", "Hypo", "Hyper", "Hypo"),
             type = c("Females", "Females", "Males", "Males", "All", "All"),
             number = c(nrow(subset(DMS_diffmeth_fem, meth.diff > 0)), nrow(subset(DMS_diffmeth_fem, meth.diff < 0)),
                        nrow(subset(DMS_diffmeth_mal, meth.diff > 0)), nrow(subset(DMS_diffmeth_mal, meth.diff < 0)),
                        nrow(subset(DMS_diffmeth_all, meth.diff > 0)), nrow(subset(DMS_diffmeth_all, meth.diff < 0))))

DMR.hyp.data <- data.frame(direction = c("Hyper", "Hypo", "Hyper", "Hypo", "Hyper", "Hypo"),
                           type = c("Females", "Females", "Males", "Males", "All", "All"),
                           number = c(nrow(subset(DMR_diffmeth_fem, meth.diff > 0)), nrow(subset(DMR_diffmeth_fem, meth.diff < 0)),
                                      nrow(subset(DMR_diffmeth_mal, meth.diff > 0)), nrow(subset(DMR_diffmeth_mal, meth.diff < 0)),
                                      nrow(subset(DMR_diffmeth_all, meth.diff > 0)), nrow(subset(DMR_diffmeth_all, meth.diff < 0))))


all.hyp.data <- data.frame(Direction = c("Hyper", "Hypo", "Hyper", "Hypo"),
                           type = c("DMS", "DMS", "DMR", "DMR"),
                               number = c(nrow(subset(DMS_diffmeth_all, meth.diff > 0)), nrow(subset(DMS_diffmeth_all, meth.diff < 0)),
                                          nrow(subset(DMR_diffmeth_all, meth.diff > 0)), nrow(subset(DMR_diffmeth_all, meth.diff < 0))))

fem.hyp.data <- data.frame(Direction = c("Hyper", "Hypo", "Hyper", "Hypo"),
                           type = c("DMS", "DMS", "DMR", "DMR"),
                           number = c(nrow(subset(DMS_diffmeth_fem, meth.diff > 0)), nrow(subset(DMS_diffmeth_fem, meth.diff < 0)),
                                      nrow(subset(DMR_diffmeth_fem, meth.diff > 0)), nrow(subset(DMR_diffmeth_fem, meth.diff < 0))))

mal.hyp.data <- data.frame(Direction = c("Hyper", "Hypo", "Hyper", "Hypo"),
                           type = c("DMS", "DMS", "DMR", "DMR"),
                           number = c(nrow(subset(DMS_diffmeth_mal, meth.diff > 0)), nrow(subset(DMS_diffmeth_mal, meth.diff < 0)),
                                      nrow(subset(DMR_diffmeth_mal, meth.diff > 0)), nrow(subset(DMR_diffmeth_mal, meth.diff < 0))))

hyp.all.plot <- plot_hyp_all(all.hyp.data)

hyp.fem.plot <- plot_hyp_all(fem.hyp.data)
hyp.mal.plot <- plot_hyp_all(mal.hyp.data)

#
# # #make panel
# hype.panel <- ggarrange(hyp.fem.plot, hyp.mal.plot, dist_fem, dist_mal,
#                         ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
# 
# hype.panel.1 <- ggarrange(hyp.fem.plot, hyp.mal.plot, 
#                         ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = TRUE, legend = "right")
# 
# hype.panel.2 <- ggarrange(dist_fem, dist_mal,
#                         ncol = 2, nrow = 1, labels = c("C", "D"), common.legend = TRUE, legend = "right")
# 
# full.panel <- ggarrange(hype.panel.1,hype.panel.2, ncol = 1, nrow = 1)
# 
# ggsave(filename = "./dev_exp/plots/finalized_tiff/panel_hyp1_dist.tiff", plot = hype.panel.1, width = 7, height = 5, units = "in", dpi = 300)
# ggsave(filename = "./dev_exp/plots/finalized_tiff/panel_hyp2_dist.tiff", plot = hype.panel.2, width = 7, height = 5, units = "in", dpi = 300)

  # #save plots 
# ggsave(filename = "./dev_exp/plots/finalized_tiff/barplots_hypo_hyper/DMS_hyp_plot.tiff", plot = hyp.DMS.plot, width = 7, height = 5, units = "in", dpi = 300)
# ggsave(filename = "./dev_exp/plots/finalized_tiff/barplots_hypo_hyper/DMR_hyp_plot.tiff", plot = hyp.DMR.plot, width = 7, height = 5, units = "in", dpi = 300)
# ggsave(filename = "./dev_exp/plots/finalized_tiff/barplots_hypo_hyper/panel_hyp_plot.tiff", plot = hype.panel, width = 8, height = 5, units = "in", dpi = 300)
ggsave(filename = "./dev_exp/plots/finalized_tiff/barplots_hypo_hyper/all_hyp_plot.tiff", plot = hyp.all.plot, width = 5, height = 5, units = "in", dpi = 300)
ggsave(filename = "./dev_exp/plots/finalized_tiff/barplots_hypo_hyper/fem_hyp_plot.tiff", plot = hyp.fem.plot, width = 5, height = 5, units = "in", dpi = 300)
ggsave(filename = "./dev_exp/plots/finalized_tiff/barplots_hypo_hyper/mal_hyp_plot.tiff", plot = hyp.mal.plot, width = 5, height = 5, units = "in", dpi = 300)

#run chi-square for sig 
chi.DMS.dat <- subset(all.hyp.data, type == "DMS")
chi.DMR.dat <- subset(all.hyp.data, type == "DMR")

fem.chi.DMS.dat <- subset(fem.hyp.data, type == "DMS")
fem.chi.DMR.dat <- subset(fem.hyp.data, type == "DMR")

mal.chi.DMS.dat <- subset(mal.hyp.data, type == "DMS")
mal.chi.DMR.dat <- subset(mal.hyp.data, type == "DMR")

chisq.test(chi.DMS.dat$number, p =c(0.5,0.5))
chisq.test(chi.DMR.dat$number, p =c(0.5,0.5))

chisq.test(fem.chi.DMS.dat$number, p =c(0.5,0.5))
chisq.test(fem.chi.DMR.dat$number, p =c(0.5,0.5))

chisq.test(mal.chi.DMS.dat$number, p =c(0.5,0.5))
chisq.test(mal.chi.DMR.dat$number, p =c(0.5,0.5))

# 
# ## Plot Manhattan plots ##
# plotManhat(DMS_diffmeth_fem, "./dev_exp/plots/finalized_tiff/man_plots/DMS_fem_manPlot.tiff")
# plotManhat(DMS_diffmeth_mal, "./dev_exp/plots/finalized_tiff/man_plots/DMS_mal_manPlot.tiff")
# plotManhat(DMS_diffmeth_all, "./dev_exp/plots/finalized_tiff/man_plots/DMS_all_manPlot.tiff")
# 
# plotManhat(DMR_diffmeth_fem, "./dev_exp/plots/finalized_tiff/man_plots/DMR_fem_manPlot.tiff")
# plotManhat(DMR_diffmeth_mal, "./dev_exp/plots/finalized_tiff/man_plots/DMR_mal_manPlot.tiff")
# plotManhat(DMR_diffmeth_all, "./dev_exp/plots/finalized_tiff/man_plots/DMR_all_manPlot.tiff")

# ## Do PCA cluster analysis ## 
#get percent methylation
DMS_fem_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMS_res/DMSdiffmeth_fem_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMS_res/DMSmeth_fem_5X.RDS")

DMS_mal_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMS_res/DMSdiffmeth_mal_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMS_res/DMSmeth_mal_5X.RDS")

DMS_all_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMS_res/DMSdiffmeth_all_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMS_res/DMSmeth_all_5X.RDS")

#get percent methylation for all CpGs
#read data files
CpG_meth_all <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSmeth_all_5X.RDS")
CpG_meth_fem <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSmeth_fem_5X.RDS")
CpG_meth_mal <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSmeth_mal_5X.RDS")

#calculate percent methylation
perc_meth_CpG_all <- as.data.frame(percMethylation(CpG_meth_all))
perc_meth_CpG_fem <- as.data.frame(percMethylation(CpG_meth_fem))
perc_meth_CpG_mal <- as.data.frame(percMethylation(CpG_meth_mal))

#add site names
CpG_meth_all <- nameSite(CpG_meth_all)
CpG_meth_fem <- nameSite(CpG_meth_fem)
CpG_meth_mal <- nameSite(CpG_meth_mal)

#add rownames to perc meth matrix
rownames(perc_meth_CpG_all) <- CpG_meth_all$site_name
rownames(perc_meth_CpG_fem) <- CpG_meth_fem$site_name
rownames(perc_meth_CpG_mal) <- CpG_meth_mal$site_name

#remove extra file
rm(CpG_meth_all)
rm(CpG_meth_fem)
rm(CpG_meth_mal)

#get perc meth for DMRs
DMR_fem_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMR_res/DMRdiffmeth_fem_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMR_res/DMR_tile_meth_fem_5X.RDS")

DMR_mal_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMR_res/DMRdiffmeth_mal_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMR_res/DMR_tile_meth_mal_5X.RDS")

DMR_all_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMR_res/DMRdiffmeth_all_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMR_res/DMR_tile_meth_all_5X.RDS")

#prep data
t_DMS_all <- prep_data(DMS_all_percMeth)
t_DMR_all <- prep_data(DMR_all_percMeth)
t_CpG_all <- prep_data(perc_meth_CpG_all)

t_DMS_fem <- prep_data(DMS_fem_percMeth)
t_DMR_fem <- prep_data(DMR_fem_percMeth)
t_CpG_fem <- prep_data(perc_meth_CpG_fem)

t_DMS_mal <- prep_data(DMS_mal_percMeth)
t_DMR_mal <- prep_data(DMR_mal_percMeth)
t_CpG_mal <- prep_data(perc_meth_CpG_mal)

#run PCA
pca_DMS_all <- prcomp(t_DMS_all, center = TRUE)
pca_DMR_all <- prcomp(t_DMR_all, center = TRUE)
pca_CpG_all <- prcomp(t_CpG_all, center = TRUE)

pca_DMS_fem <- prcomp(t_DMS_fem, center = TRUE)
pca_DMR_fem <- prcomp(t_DMR_fem, center = TRUE)
pca_CpG_fem <- prcomp(t_CpG_fem, center = TRUE)

pca_DMS_mal <- prcomp(t_DMS_mal, center = TRUE)
pca_DMR_mal <- prcomp(t_DMR_mal, center = TRUE)
pca_CpG_mal <- prcomp(t_CpG_mal, center = TRUE)

#see results
summary(pca_DMS_all)
summary(pca_DMR_all)
summary(pca_CpG_all)

summary(pca_DMS_fem)
summary(pca_DMR_fem)
summary(pca_CpG_fem)

summary(pca_DMS_mal)
summary(pca_DMR_mal)
summary(pca_CpG_mal)

#add variables to dataframe for plotting
t_DMS_all_plotDat <- add_var(t_DMS_all, pca_DMS_all)
t_DMR_all_plotDat <- add_var(t_DMR_all, pca_DMR_all)
t_CpG_all_plotDat <- add_var(t_CpG_all, pca_CpG_all)

t_DMS_fem_plotDat <- add_var(t_DMS_fem, pca_DMS_fem)
t_DMR_fem_plotDat <- add_var(t_DMR_fem, pca_DMR_fem)
t_CpG_fem_plotDat <- add_var(t_CpG_fem, pca_CpG_fem)

t_DMS_mal_plotDat <- add_var(t_DMS_mal, pca_DMS_mal)
t_DMR_mal_plotDat <- add_var(t_DMR_mal, pca_DMR_mal)
t_CpG_mal_plotDat <- add_var(t_CpG_mal, pca_CpG_mal)

#plot 
# all_DMS_pca_plot_V1 <- plot_pca(t_DMS_all_plotDat, t_DMS_all_plotDat$PC1, t_DMS_all_plotDat$PC2, "PC1 (43.48%)", "PC2 (5.75%)",
#                              "./dev_exp/plots/finalized_tiff/pca_plots/DMS_all_pca_V1.tiff")
#
# all_DMS_pca_plot_V2 <- plot_pca_V2(t_DMS_all_plotDat, t_DMS_all_plotDat$PC1, t_DMS_all_plotDat$PC2, "PC1 (43.48%)", "PC2 (5.75%)",
#                                      "./dev_exp/plots/finalized_tiff/pca_plots/DMS_all_pca_V2.tiff")
#
# all_CpG_pca_plot_V1 <- plot_pca(t_CpG_all_plotDat, t_CpG_all_plotDat$PC1, t_CpG_all_plotDat$PC2, "PC1 (11.13%)", "PC2 (0.03%)",
#                                 "./dev_exp/plots/finalized_tiff/pca_plots/CpG_all_pca_V1.tiff")
#
# all_CpG_pca_plot_V2 <- plot_pca_V2(t_CpG_all_plotDat, t_CpG_all_plotDat$PC1, t_CpG_all_plotDat$PC2, "PC1 (11.13%)", "PC2 (0.03%)",
#                                    "./dev_exp/plots/finalized_tiff/pca_plots/CpG_all_pca_V2.tiff")
#
# all_DMR_pca_plot_V1 <- plot_pca(t_DMR_all_plotDat, t_DMR_all_plotDat$PC1, t_DMR_all_plotDat$PC2, "PC1 (44.37%)", "PC2 (9.13%)",
#                                 "./dev_exp/plots/finalized_tiff/pca_plots/DMR_all_pca_V1.tiff")
#
# all_DMR_pca_plot_V2 <- plot_pca_V2(t_DMR_all_plotDat, t_DMR_all_plotDat$PC1, t_DMR_all_plotDat$PC2, "PC1 (44.37%)", "PC2 (9.13%)",
#                                    "./dev_exp/plots/finalized_tiff/pca_plots/DMR_all_pca_V2.tiff")

# fem_DMS_pca_plot <- plot_pca(t_DMS_fem_plotDat, t_DMS_fem_plotDat$PC1, t_DMS_fem_plotDat$PC2, "PC1 (42.91%)", "PC2 (8.39%)",
#                            "./dev_exp/plots/finalized_tiff/pca_plots/DMS_fem_pca.tiff")
#
# fem_CpG_pca_plot <- plot_pca(t_CpG_fem_plotDat, t_CpG_fem_plotDat$PC1, t_CpG_fem_plotDat$PC2, "PC1 (12.72%)", "PC2 (4.79%)",
#                               "./dev_exp/plots/finalized_tiff/pca_plots/CpG_fem_pca.tiff")
#
# fem_DMR_pca_plot <- plot_pca(t_DMR_fem_plotDat, t_DMR_fem_plotDat$PC1, t_DMR_fem_plotDat$PC2, "PC1 (36.83%)", "PC2 (14.42%)",
#                              "./dev_exp/plots/finalized_tiff/pca_plots/DMR_fem_pca.tiff")

# mal_DMS_pca_plot <- plot_pca(t_DMS_mal_plotDat, t_DMS_mal_plotDat$PC1, t_DMS_mal_plotDat$PC2, "PC1 (37.32%)", "PC2 (8.47%)",
#                                  "./dev_exp/plots/finalized_tiff/pca_plots/DMS_mal_pca.tiff")
#
# mal_CpG_pca_plot <- plot_pca(t_CpG_mal_plotDat, t_CpG_mal_plotDat$PC1, t_CpG_mal_plotDat$PC2, "PC1 (10.44%)", "PC2 (0.06%)",
#                              "./dev_exp/plots/finalized_tiff/pca_plots/CpG_mal_pca.tiff")
#
# mal_DMR_pca_plot <- plot_pca(t_DMR_mal_plotDat, t_DMR_mal_plotDat$PC1, t_DMR_mal_plotDat$PC2, "PC1 (40.22%)", "PC2 (6.97%)",
#                              "./dev_exp/plots/finalized_tiff/pca_plots/DMR_mal_pca.tiff")

# #make panel
# ggarrange(fem_DMS_pca_plot, fem_DMR_pca_plot, mal_DMS_pca_plot, mal_DMR_pca_plot, common.legend = TRUE, legend = "bottom",
#           labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
#
## Plot heatmaps ##
#make tank lists
AC2_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "AC2"]
AC3_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "AC3"]
AC4_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "AC4"]
AC5_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "AC5"]
AC6_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "AC6"]
AC7_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "AC7"]

C2_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "C2"]
C3_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "C3"]
C4_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "C4"]
C5_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "C5"]
C6_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "C6"]
C7_list <- t_DMS_all_plotDat$ID[t_DMS_all_plotDat$tank == "C7"]

#plot
plot_heatmap(DMS_fem_percMeth, "./dev_exp/plots/finalized_tiff/heatmaps/DMS_fem_heatmap_notanks.tiff")
plot_heatmap(DMS_mal_percMeth, "./dev_exp/plots/finalized_tiff/heatmaps/DMS_mal_heatmap_notanks.tiff")
plot_heatmap_all(DMS_all_percMeth, "./dev_exp/plots/finalized_tiff/heatmaps/DMS_all_heatmap.tiff")

plot_heatmap(DMR_fem_percMeth, "./dev_exp/plots/finalized_tiff/heatmaps/DMR_fem_heatmap_notanks.tiff")
plot_heatmap(DMR_mal_percMeth, "./dev_exp/plots/finalized_tiff/heatmaps/DMR_mal_heatmap_notanks.tiff")
plot_heatmap_all(DMR_all_percMeth, "./dev_exp/plots/finalized_tiff/heatmaps/DMR_all_heatmap.tiff")

## check for significant difference in number of sig DMSs/DMRs ##
#create contigency tables 
dms.chi.data <- data.frame("sig" = c(8769,27916), "not_sig" = c(9020131,9315923))
rownames(dms.chi.data) <- c("fem", "mal")
dms.chi.data <- t(dms.chi.data)

dmr.chi.data <- data.frame("sig" = c(51,402), "not_sig" = c(500441,501223))
rownames(dmr.chi.data) <- c("fem", "mal")
dmr.chi.data <- t(dmr.chi.data)

dms.chi <- chisq.test(dms.chi.data, correct = TRUE)
dms.chi

dmr.chi <- chisq.test(dmr.chi.data, correct = TRUE)
dmr.chi

## check for differences in variance between males and females ## 
#calculate stdev for each CpG 
fem.stdev <- apply(t_CpG_fem,2,sd)
mal.stdev <- apply(t_CpG_mal,2,sd)

#look at distributions of stdv 
hist(fem.stdev)
hist(mal.stdev)

#look at averages 
mean(fem.stdev)
mean(mal.stdev)

#do t test 
t.test(fem.stdev, mal.stdev)
