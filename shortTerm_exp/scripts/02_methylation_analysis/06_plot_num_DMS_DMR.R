#######################################################
### Goal: Plot and analyze counts of DMRs and DMSs 
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
library(UpSetR)
library(mclust)
library(lme4)

#write functions 
#function for adding site names 
nameSite <- function(data){
  data$site_name <- paste(data$chr, data$start, sep = "_")
  return(data)
}

#function for plotting barplots
barplot.funct <- function(bar.data, site, ylabel){
  ggplot(data = bar.data, aes(x=time_point, y = {{site}}, fill = sex)) + 
    geom_bar(stat="identity", position = position_dodge()) +
    theme_bw() + xlab("Time Point") + ylab(ylabel) +
    scale_fill_manual(values = c("fem" = "#481567FF", "mal" = "#9AC8CD"), name = "Sex", 
                      labels = c("all" = "All", "fem" = "Female", "mal" = "Male")) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text = element_text(size = 11),
          legend.title = element_text(size = 12)) 
}

#function for plotting upset plots
upset.plot.funct <- function(data, col){
  ggplot(data = data, aes(x = time_points)) + geom_bar(fill = col) + scale_x_upset() + theme_bw() + ylab("Intersection Size") +
    xlab("") + geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, hjust = 0.5, color = "black", size = 2) +
    theme_combmatrix(combmatrix.panel.point.color.fill = col,
                     combmatrix.panel.line.color = col,
                     combmatrix.label.text = element_text(size = 10, color = "black")) +
    theme(axis.title = element_text(size = 15, color = "black"), axis.text.y = element_text(size = 11, color = "black"))
}

#function for converting data for upset plots
convert_to_upset_format <- function(df) {
  # Extract unique time points
  unique_time_points <- unique(unlist(df$time_points))
  
  # Create an empty binary matrix
  binary_matrix <- matrix(0, nrow = nrow(df), ncol = length(unique_time_points), 
                          dimnames = list(df$site_name, unique_time_points))
  
  # Fill in binary values based on time_points
  for (i in seq_len(nrow(df))) {
    if (is.character(df$time_points[i])) {
      time_points <- as.character(df$time_points[i])
    } else {
      time_points <- unlist(df$time_points[i])
    }
    binary_matrix[i, time_points] <- 1
  }
  
  # Convert binary matrix to dataframe
  binary_df <- as.data.frame(binary_matrix)
  
  return(binary_df)
}

#function for plotting second version of upset plots
plot.upset.v2 <- function(data,col){ 
  upset(data, sets = c("0.5h", "1h", "4h", "24h", "72h"), order.by = "freq", keep.order = TRUE, 
        point.size = 2.5, line.size = 1, text.scale = c(1.5, 1.5, 1.3, 1.3, 2, 1.15), main.bar.color = col,
        matrix.color = "black", sets.bar.color = col)
}

#function to display venn diagrams
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

#create function that plots manhattan plots
plotManhat <- function(df, plotname, color_val){
  
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

#function for running chi square
chi.sq.function <- function(data){
  chisq.test(c(nrow(subset(data, meth.diff > 0)),nrow(subset(data, meth.diff < 0))), p =c(0.5,0.5))
  
}

#read data files 
DMS_diffMeth_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSdiffmeth_all_fem_5X_data.RDS")
DMS_diffMeth_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_05h/DMSdiffmeth_05h_fem_5X_data.RDS")
DMS_diffMeth_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_1h/DMSdiffmeth_1h_fem_5X_data.RDS")
DMS_diffMeth_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_4h/DMSdiffmeth_4h_fem_5X_data.RDS")
DMS_diffMeth_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_24h/DMSdiffmeth_24h_fem_5X_data.RDS")
DMS_diffMeth_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_72h/DMSdiffmeth_72h_fem_5X_data.RDS")

DMS_diffMeth_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_all/DMSdiffmeth_all_mal_5X_data.RDS")
DMS_diffMeth_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_05h/DMSdiffmeth_05h_mal_5X_data.RDS")
DMS_diffMeth_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_1h/DMSdiffmeth_1h_mal_5X_data.RDS")
DMS_diffMeth_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_4h/DMSdiffmeth_4h_mal_5X_data.RDS")
DMS_diffMeth_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_24h/DMSdiffmeth_24h_mal_5X_data.RDS")
DMS_diffMeth_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMS_res_72h/DMSdiffmeth_72h_mal_5X_data.RDS")

DMR_diffMeth_all_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_all/DMRdiffmeth_all_fem_5X_data.RDS")
DMR_diffMeth_05h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_05h/DMRdiffmeth_05h_fem_5X_data.RDS")
DMR_diffMeth_1h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_1h/DMRdiffmeth_1h_fem_5X_data.RDS")
DMR_diffMeth_4h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_4h/DMRdiffmeth_4h_fem_5X_data.RDS")
DMR_diffMeth_24h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_24h/DMRdiffmeth_24h_fem_5X_data.RDS")
DMR_diffMeth_72h_fem <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_72h/DMRdiffmeth_72h_fem_5X_data.RDS")

DMR_diffMeth_all_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_all/DMRdiffmeth_all_mal_5X_data.RDS")
DMR_diffMeth_05h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_05h/DMRdiffmeth_05h_mal_5X_data.RDS")
DMR_diffMeth_1h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_1h/DMRdiffmeth_1h_mal_5X_data.RDS")
DMR_diffMeth_4h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_4h/DMRdiffmeth_4h_mal_5X_data.RDS")
DMR_diffMeth_24h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_24h/DMRdiffmeth_24h_mal_5X_data.RDS")
DMR_diffMeth_72h_mal <- readRDS("./shortTerm_exp/data/methylkit_res/DMR_res_72h/DMRdiffmeth_72h_mal_5X_data.RDS")

#read od data files
DMS_diffMeth_all_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSdiffmeth_all_fem_5X_data.RDS")
DMS_diffMeth_05h_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_05h/DMSdiffmeth_05h_fem_5X_data.RDS")
DMS_diffMeth_1h_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_1h/DMSdiffmeth_1h_fem_5X_data.RDS")
DMS_diffMeth_4h_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_4h/DMSdiffmeth_4h_fem_5X_data.RDS")
DMS_diffMeth_24h_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_24h/DMSdiffmeth_24h_fem_5X_data.RDS")
DMS_diffMeth_72h_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_72h/DMSdiffmeth_72h_fem_5X_data.RDS")

DMS_diffMeth_all_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_all/DMSdiffmeth_all_mal_5X_data.RDS")
DMS_diffMeth_05h_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_05h/DMSdiffmeth_05h_mal_5X_data.RDS")
DMS_diffMeth_1h_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_1h/DMSdiffmeth_1h_mal_5X_data.RDS")
DMS_diffMeth_4h_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_4h/DMSdiffmeth_4h_mal_5X_data.RDS")
DMS_diffMeth_24h_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_24h/DMSdiffmeth_24h_mal_5X_data.RDS")
DMS_diffMeth_72h_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMS_res_72h/DMSdiffmeth_72h_mal_5X_data.RDS")

DMR_diffMeth_all_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMRdiffmeth_all_fem_5X_data.RDS")
DMR_diffMeth_05h_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_05h/DMRdiffmeth_05h_fem_5X_data.RDS")
DMR_diffMeth_1h_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_1h/DMRdiffmeth_1h_fem_5X_data.RDS")
DMR_diffMeth_4h_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_4h/DMRdiffmeth_4h_fem_5X_data.RDS")
DMR_diffMeth_24h_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_24h/DMRdiffmeth_24h_fem_5X_data.RDS")
DMR_diffMeth_72h_fem_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_72h/DMRdiffmeth_72h_fem_5X_data.RDS")

DMR_diffMeth_all_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_all/DMRdiffmeth_all_mal_5X_data.RDS")
DMR_diffMeth_05h_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_05h/DMRdiffmeth_05h_mal_5X_data.RDS")
DMR_diffMeth_1h_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_1h/DMRdiffmeth_1h_mal_5X_data.RDS")
DMR_diffMeth_4h_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_4h/DMRdiffmeth_4h_mal_5X_data.RDS")
DMR_diffMeth_24h_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_24h/DMRdiffmeth_24h_mal_5X_data.RDS")
DMR_diffMeth_72h_mal_od <- readRDS("./shortTerm_exp/data/methylkit_res_od/DMR_res_72h/DMRdiffmeth_72h_mal_5X_data.RDS")

#create a column with a name for each DMS/DMR
DMS_diffMeth_all_fem <- nameSite(DMS_diffMeth_all_fem)
DMS_diffMeth_05h_fem <- nameSite(DMS_diffMeth_05h_fem)
DMS_diffMeth_1h_fem <- nameSite(DMS_diffMeth_1h_fem)
DMS_diffMeth_4h_fem <- nameSite(DMS_diffMeth_4h_fem)
DMS_diffMeth_24h_fem <- nameSite(DMS_diffMeth_24h_fem)
DMS_diffMeth_72h_fem <- nameSite(DMS_diffMeth_72h_fem)

DMS_diffMeth_all_mal <- nameSite(DMS_diffMeth_all_mal)
DMS_diffMeth_05h_mal <- nameSite(DMS_diffMeth_05h_mal)
DMS_diffMeth_1h_mal <- nameSite(DMS_diffMeth_1h_mal)
DMS_diffMeth_4h_mal <- nameSite(DMS_diffMeth_4h_mal)
DMS_diffMeth_24h_mal <- nameSite(DMS_diffMeth_24h_mal)
DMS_diffMeth_72h_mal <- nameSite(DMS_diffMeth_72h_mal)

DMR_diffMeth_all_fem <- nameSite(DMR_diffMeth_all_fem)
DMR_diffMeth_05h_fem <- nameSite(DMR_diffMeth_05h_fem)
DMR_diffMeth_1h_fem <- nameSite(DMR_diffMeth_1h_fem)
DMR_diffMeth_4h_fem <- nameSite(DMR_diffMeth_4h_fem)
DMR_diffMeth_24h_fem <- nameSite(DMR_diffMeth_24h_fem)
DMR_diffMeth_72h_fem <- nameSite(DMR_diffMeth_72h_fem)

DMR_diffMeth_all_mal <- nameSite(DMR_diffMeth_all_mal)
DMR_diffMeth_05h_mal <- nameSite(DMR_diffMeth_05h_mal)
DMR_diffMeth_1h_mal <- nameSite(DMR_diffMeth_1h_mal)
DMR_diffMeth_4h_mal <- nameSite(DMR_diffMeth_4h_mal)
DMR_diffMeth_24h_mal <- nameSite(DMR_diffMeth_24h_mal)
DMR_diffMeth_72h_mal <- nameSite(DMR_diffMeth_72h_mal)

DMS_diffMeth_all_fem_od <- nameSite(DMS_diffMeth_all_fem_od)
DMS_diffMeth_05h_fem_od <- nameSite(DMS_diffMeth_05h_fem_od)
DMS_diffMeth_1h_fem_od <- nameSite(DMS_diffMeth_1h_fem_od)
DMS_diffMeth_4h_fem_od <- nameSite(DMS_diffMeth_4h_fem_od)
DMS_diffMeth_24h_fem_od <- nameSite(DMS_diffMeth_24h_fem_od)
DMS_diffMeth_72h_fem_od <- nameSite(DMS_diffMeth_72h_fem_od)

DMS_diffMeth_all_mal_od <- nameSite(DMS_diffMeth_all_mal_od)
DMS_diffMeth_05h_mal_od <- nameSite(DMS_diffMeth_05h_mal_od)
DMS_diffMeth_1h_mal_od <- nameSite(DMS_diffMeth_1h_mal_od)
DMS_diffMeth_4h_mal_od <- nameSite(DMS_diffMeth_4h_mal_od)
DMS_diffMeth_24h_mal_od <- nameSite(DMS_diffMeth_24h_mal_od)
DMS_diffMeth_72h_mal_od <- nameSite(DMS_diffMeth_72h_mal_od)

DMR_diffMeth_all_fem_od <- nameSite(DMR_diffMeth_all_fem_od)
DMR_diffMeth_05h_fem_od <- nameSite(DMR_diffMeth_05h_fem_od)
DMR_diffMeth_1h_fem_od <- nameSite(DMR_diffMeth_1h_fem_od)
DMR_diffMeth_4h_fem_od <- nameSite(DMR_diffMeth_4h_fem_od)
DMR_diffMeth_24h_fem_od <- nameSite(DMR_diffMeth_24h_fem_od)
DMR_diffMeth_72h_fem_od <- nameSite(DMR_diffMeth_72h_fem_od)

DMR_diffMeth_all_mal_od <- nameSite(DMR_diffMeth_all_mal_od)
DMR_diffMeth_05h_mal_od <- nameSite(DMR_diffMeth_05h_mal_od)
DMR_diffMeth_1h_mal_od <- nameSite(DMR_diffMeth_1h_mal_od)
DMR_diffMeth_4h_mal_od <- nameSite(DMR_diffMeth_4h_mal_od)
DMR_diffMeth_24h_mal_od <- nameSite(DMR_diffMeth_24h_mal_od)
DMR_diffMeth_72h_mal_od <- nameSite(DMR_diffMeth_72h_mal_od)

## Plot bar plots of # of DMSs and DMRs ##
#calculate number of DMSs adn DMRs
DMS.df.list <- list(DMS_diffMeth_all_fem, DMS_diffMeth_05h_fem, DMS_diffMeth_1h_fem, DMS_diffMeth_4h_fem, DMS_diffMeth_24h_fem, DMS_diffMeth_72h_fem,
               DMS_diffMeth_all_mal, DMS_diffMeth_05h_mal, DMS_diffMeth_1h_mal, DMS_diffMeth_4h_mal, DMS_diffMeth_24h_mal, DMS_diffMeth_72h_mal)

DMS <- sapply(DMS.df.list, nrow)

DMR.df.list <- list(DMR_diffMeth_all_fem, DMR_diffMeth_05h_fem, DMR_diffMeth_1h_fem, DMR_diffMeth_4h_fem, DMR_diffMeth_24h_fem, DMR_diffMeth_72h_fem,
                    DMR_diffMeth_all_mal, DMR_diffMeth_05h_mal, DMR_diffMeth_1h_mal, DMR_diffMeth_4h_mal, DMR_diffMeth_24h_mal, DMR_diffMeth_72h_mal)

DMR <- sapply(DMR.df.list, nrow)

DMS.od.df.list <- list(DMS_diffMeth_all_fem_od, DMS_diffMeth_05h_fem_od, DMS_diffMeth_1h_fem_od, DMS_diffMeth_4h_fem_od, DMS_diffMeth_24h_fem_od, DMS_diffMeth_72h_fem_od,
                    DMS_diffMeth_all_mal_od, DMS_diffMeth_05h_mal_od, DMS_diffMeth_1h_mal_od, DMS_diffMeth_4h_mal_od, DMS_diffMeth_24h_mal_od, DMS_diffMeth_72h_mal_od)

DMS.od <- sapply(DMS.od.df.list, nrow)

DMR.od.df.list <- list(DMR_diffMeth_all_fem_od, DMR_diffMeth_05h_fem_od, DMR_diffMeth_1h_fem_od, DMR_diffMeth_4h_fem_od, DMR_diffMeth_24h_fem_od, DMR_diffMeth_72h_fem_od,
                    DMR_diffMeth_all_mal_od, DMR_diffMeth_05h_mal_od, DMR_diffMeth_1h_mal_od, DMR_diffMeth_4h_mal_od, DMR_diffMeth_24h_mal_od, DMR_diffMeth_72h_mal_od)

DMR.od <- sapply(DMR.od.df.list, nrow)

#make data frame with number of sites
barplot.data <- data.frame(num_DMS = DMS, num_DMR = DMR, 
                           time_point = c("all", "0.5h", "1h", "4h", "24h", "72h",
                                          "all", "0.5h", "1h", "4h", "24h", "72h"),
                           sex = c("fem", "fem", "fem", "fem", "fem", "fem",
                                   "mal", "mal", "mal", "mal", "mal", "mal"))

barplot.data$time_point <- factor(barplot.data$time_point, ordered = TRUE, 
                                  levels = c("0.5h", "1h", "4h", "24h", "72h", "all"))

barplot.od.data <- data.frame(num_DMS = DMS.od, num_DMR = DMR.od, 
                           time_point = c("all", "0.5h", "1h", "4h", "24h", "72h",
                                          "all", "0.5h", "1h", "4h", "24h", "72h"),
                           sex = c("fem", "fem", "fem", "fem", "fem", "fem",
                                   "mal", "mal", "mal", "mal", "mal", "mal"))

barplot.od.data$time_point <- factor(barplot.od.data$time_point, ordered = TRUE, 
                                  levels = c("0.5h", "1h", "4h", "24h", "72h", "all"))

#plot barplots 
DMS.barplot <- barplot.funct(barplot.data, num_DMS, "Number of DMSs")
DMR.barplot <- barplot.funct(barplot.data, num_DMR, "Number of DMRs")

DMS.od.barplot <- barplot.funct(barplot.od.data, num_DMS, "Number of DMSs")
DMR.od.barplot <- barplot.funct(barplot.od.data, num_DMR, "Number of DMRs")

#make a panel 
barplot.panel <- ggarrange(DMS.barplot, DMR.barplot, labels = c("A", "B"),
                           nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom") 

barplot.od.panel <- ggarrange(DMS.od.barplot, DMR.od.barplot, labels = c("A", "B"),
                           nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom") 

#save plots
tiff("./shortTerm_exp/plots/finalized_tiff/DMSbarPlot.tiff", units="in", width = 6, height = 4, res = 600)
DMS.barplot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMRbarPlot.tiff", units="in", width = 6, height = 4, res = 600)
DMR.barplot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/panel_barPlot.tiff", units="in", width = 6, height = 4, res = 600)
barplot.panel
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMSbarPlot_od.tiff", units="in", width = 6, height = 4, res = 600)
DMS.od.barplot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMRbarPlot_od.tiff", units="in", width = 6, height = 4, res = 600)
DMR.od.barplot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/panel_od_barPlot.tiff", units="in", width = 6, height = 4, res = 600)
barplot.od.panel
dev.off()

## Plot upset plots ## 
#add time points 
DMS_diffMeth_05h_fem$time_point <- "0.5h"
DMS_diffMeth_1h_fem$time_point <- "1h"
DMS_diffMeth_4h_fem$time_point <- "4h"
DMS_diffMeth_24h_fem$time_point <- "24h"
DMS_diffMeth_72h_fem$time_point <- "72h"

DMS_diffMeth_05h_mal$time_point <- "0.5h"
DMS_diffMeth_1h_mal$time_point <- "1h"
DMS_diffMeth_4h_mal$time_point <- "4h"
DMS_diffMeth_24h_mal$time_point <- "24h"
DMS_diffMeth_72h_mal$time_point <- "72h"

DMR_diffMeth_05h_fem$time_point <- "0.5h"
DMR_diffMeth_1h_fem$time_point <- "1h"
DMR_diffMeth_4h_fem$time_point <- "4h"
DMR_diffMeth_24h_fem$time_point <- "24h"
DMR_diffMeth_72h_fem$time_point <- "72h"

DMR_diffMeth_05h_mal$time_point <- "0.5h"
DMR_diffMeth_1h_mal$time_point <- "1h"
DMR_diffMeth_4h_mal$time_point <- "4h"
DMR_diffMeth_24h_mal$time_point <- "24h"
DMR_diffMeth_72h_mal$time_point <- "72h"

#format data frames for upset plots 
DMS.fem.upset <- rbind(DMS_diffMeth_05h_fem,DMS_diffMeth_1h_fem,DMS_diffMeth_4h_fem,DMS_diffMeth_24h_fem,DMS_diffMeth_72h_fem)
DMS.fem.upset <- DMS.fem.upset %>% group_by(site_name) %>% summarize(time_points = list(time_point))

DMS.mal.upset <- rbind(DMS_diffMeth_05h_mal,DMS_diffMeth_1h_mal,DMS_diffMeth_4h_mal,DMS_diffMeth_24h_mal,DMS_diffMeth_72h_mal)
DMS.mal.upset <- DMS.mal.upset %>% group_by(site_name) %>% summarize(time_points = list(time_point))

DMR.fem.upset <- rbind(DMR_diffMeth_05h_fem,DMR_diffMeth_1h_fem,DMR_diffMeth_4h_fem,DMR_diffMeth_24h_fem,DMR_diffMeth_72h_fem)
DMR.fem.upset <- DMR.fem.upset %>% group_by(site_name) %>% summarize(time_points = list(time_point))

DMR.mal.upset <- rbind(DMR_diffMeth_05h_mal,DMR_diffMeth_1h_mal,DMR_diffMeth_4h_mal,DMR_diffMeth_24h_mal,DMR_diffMeth_72h_mal)
DMR.mal.upset <- DMR.mal.upset %>% group_by(site_name) %>% summarize(time_points = list(time_point))

#plot upsets
DMS.fem.upset.plot <- upset.plot.funct(DMS.fem.upset,"#481567FF")
DMS.mal.upset.plot <- upset.plot.funct(DMS.mal.upset, "#9AC8CD")

DMR.fem.upset.plot <- upset.plot.funct(DMR.fem.upset,"#481567FF")
DMR.mal.upset.plot <- upset.plot.funct(DMR.mal.upset,"#9AC8CD")

#save plots
tiff("./shortTerm_exp/plots/finalized_tiff/DMS_fem_upset.tiff", units="in", width = 7, height = 6, res = 600)
DMS.fem.upset.plot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMS_mal_upset.tiff", units="in", width = 7, height = 6, res = 600)
DMS.mal.upset.plot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMR_fem_upset.tiff", units="in", width = 7, height = 6, res = 600)
DMR.fem.upset.plot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMR_mal_upset.tiff", units="in", width = 7, height = 6, res = 600)
DMR.mal.upset.plot
dev.off()

#try different way of plotting 
#convert data 
DMS.fem.upset.new <- convert_to_upset_format(DMS.fem.upset)
DMS.mal.upset.new <- convert_to_upset_format(DMS.mal.upset)

DMR.fem.upset.new <- convert_to_upset_format(DMR.fem.upset)
DMR.mal.upset.new <- convert_to_upset_format(DMR.mal.upset)

#plot upsets 
DMS.fem.upset.new.plot <- plot.upset.v2(DMS.fem.upset.new, "#481567FF")
DMS.mal.upset.new.plot <- plot.upset.v2(DMS.mal.upset.new, "#9AC8CD")

DMR.fem.upset.new.plot <- plot.upset.v2(DMR.fem.upset.new, "#481567FF")
DMR.mal.upset.new.plot <- plot.upset.v2(DMR.mal.upset.new, "#9AC8CD")

#save plots
tiff("./shortTerm_exp/plots/finalized_tiff/DMS_fem_upset_v2.tiff", units="in", width = 10, height = 6, res = 600)
DMS.fem.upset.new.plot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMS_mal_upset_v2.tiff", units="in", width = 10, height = 6, res = 600)
DMS.mal.upset.new.plot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMS_all_upset_v2.tiff", units="in", width = 10, height = 6, res = 600)
DMS.all.upset.new.plot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMR_fem_upset_v2.tiff", units="in", width = 10, height = 6, res = 600)
DMR.fem.upset.new.plot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMR_mal_upset_v2.tiff", units="in", width = 10, height = 6, res = 600)
DMR.mal.upset.new.plot
dev.off()

tiff("./shortTerm_exp/plots/finalized_tiff/DMR_all_upset_v2.tiff", units="in", width = 10, height = 6, res = 600)
DMR.all.upset.new.plot
dev.off()

#find overlapping DMRs in all time points
Reduce(intersect, list(DMS_diffMeth_05h_fem$site_name,DMS_diffMeth_1h_fem$site_name,DMS_diffMeth_4h_fem$site_name,
                       DMS_diffMeth_24h_fem$site_name,DMS_diffMeth_72h_fem$site_name))

Reduce(intersect, list(DMS_diffMeth_05h_mal$site_name,DMS_diffMeth_1h_mal$site_name,DMS_diffMeth_4h_mal$site_name,
                       DMS_diffMeth_24h_mal$site_name,DMS_diffMeth_72h_mal$site_name))

Reduce(intersect, list(DMR_diffMeth_05h_fem$site_name,DMR_diffMeth_1h_fem$site_name,DMR_diffMeth_4h_fem$site_name,
                       DMR_diffMeth_24h_fem$site_name,DMR_diffMeth_72h_fem$site_name))

Reduce(intersect, list(DMR_diffMeth_05h_mal$site_name,DMR_diffMeth_1h_mal$site_name,DMR_diffMeth_4h_mal$site_name,
                       DMR_diffMeth_24h_mal$site_name,DMR_diffMeth_72h_mal$site_name))

#make panels
upset.panel <- ggarrange(DMS.fem.upset.plot, DMS.mal.upset.plot, DMR.fem.upset.plot, DMR.mal.upset.plot,
                         labels = c("A", "B", "C", "D"), nrow = 2, ncol = 2)

tiff("./shortTerm_exp/plots/finalized_tiff/upset_panel.tiff", units="in", width = 13, height = 8, res = 600)
upset.panel
dev.off()

## Plot venn diagrams ##
# Plot one for males vs females across all time points
#DMS
DMS.mal <- DMS_diffMeth_all_mal$site_name
DMS.fem <- DMS_diffMeth_all_fem$site_name

DMS.femvmal <- list(DMS.mal, DMS.fem)

venn_object <- venn.diagram(DMS.femvmal, category.names = c("Males", "Females"),
             filename = "./shortTerm_exp/plots/finalized_tiff/venn_diagrams/DMS_femVmal_venn_od.png",
             imagetype = "tiff", height = 2000, width = 2000, resolution = 600,
             lwd = 2, col = c("#1F968BFF", "#481567FF"),
             fill = c(alpha("#1F968BFF", 0.8), alpha("#481567FF", 0.8)),
             cex = 1, cat.cex = 1, scaled = FALSE, fontfamily = "sans",
             cat.fontfamily = "sans",cat.default.pos = "outer",cat.pos = c(-27, 27),
             cat.dist = c(0.055, 0.055))

#DMRs
DMR.mal <- DMR_diffMeth_all_mal$site_name
DMR.fem <- DMR_diffMeth_all_fem$site_name

DMR.femvmal <- list(DMR.mal, DMR.fem)

venn_object <- venn.diagram(DMR.femvmal, category.names = c("Males", "Females"),
                            filename = "./shortTerm_exp/plots/finalized_tiff/venn_diagrams/DMR_femVmal_venn_od.png",
                            imagetype = "tiff", height = 2000, width = 2000, resolution = 600,
                            lwd = 2, col = c("#1F968BFF", "#481567FF"),
                            fill = c(alpha("#1F968BFF", 0.8), alpha("#481567FF", 0.8)),
                            cex = 1, cat.cex = 1, scaled = FALSE, fontfamily = "sans",
                            cat.fontfamily = "sans",cat.default.pos = "outer",cat.pos = c(-27, 27),
                            cat.dist = c(0.055, 0.055), inverted = TRUE) #note have to invert so that it is in the same order as the previous one

## Plot Manhattan plots ##
color_val1 <- rep(c("#1F968BFF", "#39568CFF"), length.out = 23)

DMS.man.all.fem <- plotManhat(DMS_diffMeth_all_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_allTP_fem_manPlot_od.tiff", color_val1)
DMS.man.fem.05h <- plotManhat(DMS_diffMeth_05h_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_05h_fem_manPlot_od.tiff", color_val1)
DMS.man.fem.1h <- plotManhat(DMS_diffMeth_1h_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_1h_fem_manPlot_od.tiff", color_val1)
DMS.man.fem.4h <- plotManhat(DMS_diffMeth_4h_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_4h_fem_manPlot_od.tiff", color_val1)
DMS.man.fem.24h <- plotManhat(DMS_diffMeth_24h_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_24h_fem_manPlot_od.tiff", color_val1)
DMS.man.fem.72h <- plotManhat(DMS_diffMeth_72h_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_72h_fem_manPlot_od.tiff", color_val1)

DMS.man.all.mal <- plotManhat(DMS_diffMeth_all_mal, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_allTP_mal_manPlot_od.tiff", color_val1)
DMS.man.mal.05h <- plotManhat(DMS_diffMeth_05h_mal, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_05h_mal_manPlot_od.tiff", color_val1)
DMS.man.mal.1h <- plotManhat(DMS_diffMeth_1h_mal, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_1h_mal_manPlot_od.tiff", color_val1)
DMS.man.mal.4h <- plotManhat(DMS_diffMeth_4h_mal, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_4h_mal_manPlot_od.tiff", color_val1)
DMS.man.mal.24h <- plotManhat(DMS_diffMeth_24h_mal, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_24h_mal_manPlot_od.tiff", color_val1)
DMS.man.mal.72h <- plotManhat(DMS_diffMeth_72h_mal, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_72h_mal_manPlot_od.tiff", color_val1)

DMR.man.all.fem <- plotManhat(DMR_diffMeth_all, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_allTP_allSex_manPlot_od.tiff", color_val1)
DMR.man.fem.05h <- plotManhat(DMR_diffMeth_05h, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_05h_allSex_manPlot_od.tiff", color_val1)
DMR.man.fem.1h <- plotManhat(DMR_diffMeth_1h, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_1h_allSex_manPlot_od.tiff", color_val1)
DMR.man.fem.4h <- plotManhat(DMR_diffMeth_4h, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_4h_allSex_manPlot_od.tiff", color_val1)
DMR.man.fem.24h <- plotManhat(DMR_diffMeth_24h, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_24h_allSex_manPlot_od.tiff", color_val1)
DMR.man.fem.72h <- plotManhat(DMR_diffMeth_72h, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_72h_allSex_manPlot_od.tiff", color_val1)

DMR.man.all.mal <- plotManhat(DMR_diffMeth_all_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_allTP_fem_manPlot_od.tiff", color_val1)
DMR.man.mal.05h <- plotManhat(DMR_diffMeth_05h_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_05h_fem_manPlot_od.tiff", color_val1)
DMR.man.mal.1h <- plotManhat(DMR_diffMeth_1h_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_1h_fem_manPlot_od.tiff", color_val1)
DMR.man.mal.4h <- plotManhat(DMR_diffMeth_4h_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_4h_fem_manPlot_od.tiff", color_val1)
DMR.man.mal.24h <- plotManhat(DMR_diffMeth_24h_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_24h_fem_manPlot_od.tiff", color_val1)
DMR.man.mal.72h <- plotManhat(DMR_diffMeth_72h_fem, "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_72h_fem_manPlot_od.tiff", color_val1)

#make panels
DMS.man.panel.fem <- ggarrange(DMS.man.fem.05h, DMS.man.fem.1h, DMS.man.fem.4h, DMS.man.fem.24h, DMS.man.fem.72h, labels = c("A", "B", "C", "D", "E"),
                               nrow = 3, ncol = 2)
ggsave(filename = "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_fem_man_panel_od.tiff", plot = DMS.man.panel.fem, width = 12, height = 12, units = "in", dpi = 300)

DMS.man.panel.mal <- ggarrange(DMS.man.mal.05h, DMS.man.mal.1h, DMS.man.mal.4h, DMS.man.mal.24h, DMS.man.mal.72h, labels = c("A", "B", "C", "D", "E"),
                               nrow = 3, ncol = 2)
ggsave(filename = "./shortTerm_exp/plots/finalized_tiff/man_plots/DMS_mal_man_panel_od.tiff", plot = DMS.man.panel.mal, width = 12, height = 12, units = "in", dpi = 300)


DMR.man.panel.fem <- ggarrange(DMR.man.fem.05h, DMR.man.fem.1h, DMR.man.fem.4h, DMR.man.fem.24h, DMR.man.fem.72h, labels = c("A", "B", "C", "D", "E"),
                                nrow = 3, ncol = 2)
ggsave(filename = "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_fem_man_panel_od.tiff", plot = DMR.man.panel.fem, width = 12, height = 12, units = "in", dpi = 300)

DMR.man.panel.mal <- ggarrange(DMR.man.mal.05h, DMR.man.mal.1h, DMR.man.mal.4h, DMR.man.mal.24h, DMR.man.mal.72h, labels = c("A", "B", "C", "D", "E"),
                               nrow = 3, ncol = 2)
ggsave(filename = "./shortTerm_exp/plots/finalized_tiff/man_plots/DMR_mal_man_panel_od.tiff", plot = DMR.man.panel.mal, width = 12, height = 12, units = "in", dpi = 300)

## Assess changes in methylation ##
#check max and min methylation changes in DMRs
min(abs(DMR_diffMeth_all_fem$meth.diff))
max(abs(DMR_diffMeth_all_fem$meth.diff))

min(abs(DMR_diffMeth_all_mal$meth.diff))
max(abs(DMR_diffMeth_all_mal$meth.diff))

min(abs(rbind(DMR_diffMeth_05h_fem,DMR_diffMeth_1h_fem,DMR_diffMeth_4h_fem,DMR_diffMeth_24h_fem,DMR_diffMeth_72h_fem)$meth.diff))
max(abs(rbind(DMR_diffMeth_05h_fem,DMR_diffMeth_1h_fem,DMR_diffMeth_4h_fem,DMR_diffMeth_24h_fem,DMR_diffMeth_72h_fem)$meth.diff))

min(abs(rbind(DMR_diffMeth_05h_mal,DMR_diffMeth_1h_mal,DMR_diffMeth_4h_mal,DMR_diffMeth_24h_mal,DMR_diffMeth_72h_mal)$meth.diff))
max(abs(rbind(DMR_diffMeth_05h_mal,DMR_diffMeth_1h_mal,DMR_diffMeth_4h_mal,DMR_diffMeth_24h_mal,DMR_diffMeth_72h_mal)$meth.diff))

#check hyper vs hypo methylation
#pooled time points
nrow(subset(DMS_diffMeth_all_fem, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_all_fem, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_all_fem, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_all_fem, meth.diff > 0)) #hyper

nrow(subset(DMS_diffMeth_all_mal, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_all_mal, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_all_mal, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_all_mal, meth.diff > 0)) #hyper

#female unpooled time points
nrow(subset(DMS_diffMeth_05h_fem, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_05h_fem, meth.diff > 0)) #hyper

nrow(subset(DMS_diffMeth_1h_fem, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_1h_fem, meth.diff > 0)) #hyper

nrow(subset(DMS_diffMeth_4h_fem, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_4h_fem, meth.diff > 0)) #hyper

nrow(subset(DMS_diffMeth_24h_fem, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_24h_fem, meth.diff > 0)) #hyper

nrow(subset(DMS_diffMeth_72h_fem, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_72h_fem, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_05h_fem, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_05h_fem, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_1h_fem, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_1h_fem, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_4h_fem, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_4h_fem, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_24h_fem, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_24h_fem, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_72h_fem, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_72h_fem, meth.diff > 0)) #hyper

#male unpooled time points
nrow(subset(DMS_diffMeth_05h_mal, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_05h_mal, meth.diff > 0)) #hyper

nrow(subset(DMS_diffMeth_1h_mal, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_1h_mal, meth.diff > 0)) #hyper

nrow(subset(DMS_diffMeth_4h_mal, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_4h_mal, meth.diff > 0)) #hyper

nrow(subset(DMS_diffMeth_24h_mal, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_24h_mal, meth.diff > 0)) #hyper

nrow(subset(DMS_diffMeth_72h_mal, meth.diff < 0)) #hypo
nrow(subset(DMS_diffMeth_72h_mal, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_05h_mal, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_05h_mal, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_1h_mal, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_1h_mal, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_4h_mal, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_4h_mal, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_24h_mal, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_24h_mal, meth.diff > 0)) #hyper

nrow(subset(DMR_diffMeth_72h_mal, meth.diff < 0)) #hypo
nrow(subset(DMR_diffMeth_72h_mal, meth.diff > 0)) #hyper

#run chi square tests 
#DMS
chisq.test(c(nrow(subset(DMS_diffMeth_all_fem, meth.diff > 0)),nrow(subset(DMS_diffMeth_all_fem, meth.diff < 0))), p =c(0.5,0.5))
chisq.test(c(nrow(subset(DMS_diffMeth_all_mal, meth.diff > 0)),nrow(subset(DMS_diffMeth_all_mal, meth.diff < 0))), p =c(0.5,0.5))

chi.sq.function(DMS_diffMeth_05h_fem)
chi.sq.function(DMS_diffMeth_1h_fem)
chi.sq.function(DMS_diffMeth_4h_fem)
chi.sq.function(DMS_diffMeth_24h_fem)
chi.sq.function(DMS_diffMeth_72h_fem)

chi.sq.function(DMR_diffMeth_05h_fem)
chi.sq.function(DMR_diffMeth_1h_fem)
chi.sq.function(DMR_diffMeth_4h_fem)
chi.sq.function(DMR_diffMeth_24h_fem)
chi.sq.function(DMR_diffMeth_72h_fem)

chi.sq.function(DMS_diffMeth_05h_mal)
chi.sq.function(DMS_diffMeth_1h_mal)
chi.sq.function(DMS_diffMeth_4h_mal)
chi.sq.function(DMS_diffMeth_24h_mal)
chi.sq.function(DMS_diffMeth_72h_mal)

chi.sq.function(DMR_diffMeth_05h_mal)
chi.sq.function(DMR_diffMeth_1h_mal)
chi.sq.function(DMR_diffMeth_4h_mal)
chi.sq.function(DMR_diffMeth_24h_mal)
chi.sq.function(DMR_diffMeth_72h_mal)

#DMR
chisq.test(c(nrow(subset(DMR_diffMeth_all_fem, meth.diff > 0)),nrow(subset(DMR_diffMeth_all_fem, meth.diff < 0))), p =c(0.5,0.5))
chisq.test(c(nrow(subset(DMR_diffMeth_all_mal, meth.diff > 0)),nrow(subset(DMR_diffMeth_all_mal, meth.diff < 0))), p =c(0.5,0.5))

#save lists of DMRs and DMSs for use in future scripts
saveRDS(DMS_diffMeth_05h_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMS_05h_fem_od.RDS")
saveRDS(DMS_diffMeth_1h_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMS_1h_fem_od.RDS")
saveRDS(DMS_diffMeth_4h_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMS_4h_fem_od.RDS")
saveRDS(DMS_diffMeth_24h_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMS_24h_fem_od.RDS")
saveRDS(DMS_diffMeth_72h_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMS_72h_fem_od.RDS")
saveRDS(DMS_diffMeth_all_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMS_allTimePoints_fem_od.RDS")

saveRDS(DMS_diffMeth_05h_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMS_05h_mal_od.RDS")
saveRDS(DMS_diffMeth_1h_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMS_1h_mal_od.RDS")
saveRDS(DMS_diffMeth_4h_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMS_4h_mal_od.RDS")
saveRDS(DMS_diffMeth_24h_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMS_24h_mal_od.RDS")
saveRDS(DMS_diffMeth_72h_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMS_72h_mal_od.RDS")
saveRDS(DMS_diffMeth_all_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMS_allTimePoints_mal_od.RDS")

saveRDS(DMR_diffMeth_05h_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMR_05h_fem_od.RDS")
saveRDS(DMR_diffMeth_1h_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMR_1h_fem_od.RDS")
saveRDS(DMR_diffMeth_4h_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMR_4h_fem_od.RDS")
saveRDS(DMR_diffMeth_24h_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMR_24h_fem_od.RDS")
saveRDS(DMR_diffMeth_72h_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMR_72h_fem_od.RDS")
saveRDS(DMR_diffMeth_all_fem$site_name, file = "./shortTerm_exp/data/site_lists/DMR_allTimePoints_fem_od.RDS")

saveRDS(DMR_diffMeth_05h_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMR_05h_mal_od.RDS")
saveRDS(DMR_diffMeth_1h_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMR_1h_mal_od.RDS")
saveRDS(DMR_diffMeth_4h_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMR_4h_mal_od.RDS")
saveRDS(DMR_diffMeth_24h_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMR_24h_mal_od.RDS")
saveRDS(DMR_diffMeth_72h_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMR_72h_mal_od.RDS")
saveRDS(DMR_diffMeth_all_mal$site_name, file = "./shortTerm_exp/data/site_lists/DMR_allTimePoints_mal_od.RDS")
