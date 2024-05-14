#######################################################
### Goal: Plot and analyze sex diff results data
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
library(methylKit)
library(genomation)
library(DescTools)
library(tibble)

# load data 
diffmethChr.dev <- readRDS("./sexDiff/methylKit_res/onlyControl/DMSdiffmethchr5x_dev.RDS")
diffmethChr.st <- readRDS("./sexDiff/methylKit_res/onlyControl/DMSdiffmethchr5x_st.RDS")
diffmethChr.st.all <- readRDS("./sexDiff/methylKit_res/DMSdiffmethchr5x_st_all.RDS")

DMS.diffmeth.dev <- readRDS("./sexDiff/methylKit_res/onlyControl/DMSdiffmeth5x_dev_data.RDS")
DMS.diffmeth.st <- readRDS("./sexDiff/methylKit_res/onlyControl/DMSdiffmeth5x_st_data.RDS")
DMS.diffmeth.st.all <- readRDS("./sexDiff/methylKit_res/DMSdiffmeth5x_st_all_data.RDS")

DMR.diffmeth.dev <- readRDS("./sexDiff/methylKit_res/onlyControl/DMRdiffmeth5x_data_dev.RDS")
DMR.diffmeth.st <- readRDS("./sexDiff/methylKit_res/onlyControl/DMRdiffmeth5x_data_st.RDS")
DMR.diffmeth.st.all <- readRDS("./sexDiff/methylKit_res/DMRdiffmeth5x_data_st_all.RDS")

DMS.anno.dev <- readRDS("./sexDiff/methylKit_res/onlyControl/DMS_annStats_num_dev.RDS")
DMS.anno.st <- readRDS("./sexDiff/methylKit_res/onlyControl/DMS_annStats_num_st.RDS")
DMS.anno.st.all <- readRDS("./sexDiff/methylKit_res/DMS_annStats_num_st_all.RDS")

DMR.anno.dev <- readRDS("./sexDiff/methylKit_res/onlyControl/DMR_annStats_num_dev.RDS")
DMR.anno.st <- readRDS("./sexDiff/methylKit_res/onlyControl/DMR_annStats_num_st.RDS")
DMR.anno.st.all <- readRDS("./sexDiff/methylKit_res/DMR_annStats_num_st_all.RDS")

cpg.anno.dev <- readRDS("./sexDiff/methylKit_res/onlyControl/CpG_annStats_num_dev.RDS")
cpg.anno.st <- readRDS("./sexDiff/methylKit_res/onlyControl/CpG_annStats_num_st.RDS")
cpg.anno.st.all <- readRDS("./sexDiff/methylKit_res/CpG_annStats_num_st_all.RDS")

DMS.tss.dev <- readRDS("./sexDiff/methylKit_res/onlyControl/DMS_TSS_dev.RDS")
DMS.tss.st <- readRDS("./sexDiff/methylKit_res/onlyControl/DMS_TSS_st.RDS")
DMS.tss.st.all <- readRDS("./sexDiff/methylKit_res/DMS_TSS_st_all.RDS")

## Analyse and plot Manhattans ##
#calculate percentage of DMS on sex chromosome
table(DMS.diffmeth.dev$chr=="NC_024342.1")
396/(2357+396) * 100

table(DMR.diffmeth.dev$chr=="NC_024342.1")
1/(1+17) * 100

table(DMS.diffmeth.st$chr=="NC_024342.1")
1228/(3728+1228) * 100
 
table(DMR.diffmeth.st$chr=="NC_024342.1")
13/(22+13) * 100

table(DMS.diffmeth.st.all$chr=="NC_024342.1")
344/(344+185) * 100

table(DMR.diffmeth.st.all$chr=="NC_024342.1")
#100%

#rename chromosomes
rename.chr <- function(df){
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
              return(df)
}

DMS.diffmeth.dev <- rename.chr(DMS.diffmeth.dev)
DMS.diffmeth.st <- rename.chr(DMS.diffmeth.st)
DMS.diffmeth.st.all <- rename.chr(DMS.diffmeth.st.all)

DMR.diffmeth.dev$chr <- as.character(DMR.diffmeth.dev$chr)
DMR.diffmeth.st$chr <- as.character(DMR.diffmeth.st$chr)
DMR.diffmeth.st.all$chr <- as.character(DMR.diffmeth.st.all$chr)

DMR.diffmeth.dev <- rename.chr(DMR.diffmeth.dev)
DMR.diffmeth.st <- rename.chr(DMR.diffmeth.st)
DMR.diffmeth.st.all <- rename.chr(DMR.diffmeth.st.all)

#convert to numeric
DMS.diffmeth.dev$chr <- as.numeric(DMS.diffmeth.dev$chr)
DMS.diffmeth.st$chr <- as.numeric(DMS.diffmeth.st$chr)
DMS.diffmeth.st.all$chr <- as.numeric(DMS.diffmeth.st.all$chr)

DMR.diffmeth.dev$chr <- as.numeric(DMR.diffmeth.dev$chr)
DMR.diffmeth.st$chr <- as.numeric(DMR.diffmeth.st$chr)
DMR.diffmeth.st.all$chr <- as.numeric(DMR.diffmeth.st.all$chr)

#add a column with a name for each DMS 
DMS.diffmeth.dev$DMS <- paste(DMS.diffmeth.dev$chr, DMS.diffmeth.dev$start, sep = "_")
DMS.diffmeth.st$DMS <- paste(DMS.diffmeth.st$chr, DMS.diffmeth.st$start, sep = "_")
DMS.diffmeth.st.all$DMS <- paste(DMS.diffmeth.st.all$chr, DMS.diffmeth.st.all$start, sep = "_")

DMR.diffmeth.dev$DMR <- paste(DMR.diffmeth.dev$chr, DMR.diffmeth.dev$start, sep = "_")
DMR.diffmeth.st$DMR <- paste(DMR.diffmeth.st$chr, DMR.diffmeth.st$start, sep = "_")
DMR.diffmeth.st.all$DMR <- paste(DMR.diffmeth.st.all$chr, DMR.diffmeth.st.all$start, sep = "_")

#compute cumulative position of SNP
DMS.dev.don <- DMS.diffmeth.dev %>%
  #compute chr size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>%
  #calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  #add this info to the initial dataset
  left_join(DMS.diffmeth.dev, ., by=c("chr"="chr")) %>%
  #add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate(DMScum=start+tot)

DMS.st.don <- DMS.diffmeth.st %>%
  #compute chr size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>%
  #calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  #add this info to the initial dataset
  left_join(DMS.diffmeth.st, ., by=c("chr"="chr")) %>%
  #add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate(DMScum=start+tot)

DMS.st.all.don <- DMS.diffmeth.st.all %>%
  #compute chr size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>%
  #calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  #add this info to the initial dataset
  left_join(DMS.diffmeth.st.all, ., by=c("chr"="chr")) %>%
  #add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate(DMScum=start+tot)

DMR.dev.don <- DMR.diffmeth.dev %>%
  #compute chr size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>%
  #calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  #add this info to the initial dataset
  left_join(DMR.diffmeth.dev, ., by=c("chr"="chr")) %>%
  #add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate(DMScum=start+tot)

DMR.st.don <- DMR.diffmeth.st %>%
  #compute chr size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>%
  #calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  #add this info to the initial dataset
  left_join(DMR.diffmeth.st, ., by=c("chr"="chr")) %>%
  #add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate(DMScum=start+tot)

DMR.st.all.don <- DMR.diffmeth.st.all %>%
  #compute chr size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>%
  #calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  #add this info to the initial dataset
  left_join(DMR.diffmeth.st.all, ., by=c("chr"="chr")) %>%
  #add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate(DMScum=start+tot)

#create axis
DMS.axisDev <- DMS.dev.don %>% group_by(chr) %>% summarize(center=(max(DMScum) + min(DMScum))/2)
DMS.axisSt <- DMS.st.don %>% group_by(chr) %>% summarize(center=(max(DMScum) + min(DMScum))/2)
DMS.axisSt.all <- DMS.st.all.don %>% group_by(chr) %>% summarize(center=(max(DMScum) + min(DMScum))/2)

DMR.axisDev <- DMR.dev.don %>% group_by(chr) %>% summarize(center=(max(DMScum) + min(DMScum))/2)
DMR.axisSt <- DMR.st.don %>% group_by(chr) %>% summarize(center=(max(DMScum) + min(DMScum))/2)
DMR.axisSt.all <- DMR.st.all.don %>% group_by(chr) %>% summarize(center=(max(DMScum) + min(DMScum))/2)

#create lists of colors for chromosomes
color_val_dev <- rep(c("grey", "skyblue"), length.out = 23)
color_val_dev[12] <- "orange"

color_val_st <- rep(c("#1F968BFF", "#39568CFF"), length.out = 23)
color_val_st[12] <- "#B8DE29FF"

#make manhattan plots 
plot.manhat <- function(data, axisDF, color_val){
ggplot(data, aes(x=DMScum, y=meth.diff)) +
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = color_val) +
  scale_x_continuous(label = axisDF$chr, breaks = axisDF$center) +
  theme_light() + ylab("Change in % methylation") + xlab("Chromosome") +
  geom_hline(yintercept = 0, color = "red", linetype="dashed") +
  theme( 
    legend.position="none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
}

# DMS.dev.man <- plot.manhat(DMS.dev.don, DMS.axisDev, color_val_dev)
# DMS.dev.man 
# 
# tiff("./sexDiff/plots/DMS_sexDevManPlot.tiff", units="in", width = 7, height = 4, res = 600)
# DMS.dev.man
# dev.off()
# 
DMS.st.man <- plot.manhat(DMS.st.don, DMS.axisSt, color_val_st)
DMS.st.man

tiff("./sexDiff/plots/DMS_sexStManPlot.tiff", units="in", width = 7, height = 4, res = 600)
DMS.st.man
dev.off()
# 
# DMR.dev.man <- plot.manhat(DMR.dev.don, DMR.axisDev, color_val_dev)
# DMR.dev.man 
# 
# tiff("./sexDiff/plots/DMR_sexDevManPlot.tiff", units="in", width = 7, height = 4, res = 600)
# DMR.dev.man
# dev.off()
# 
# DMR.st.man <- plot.manhat(DMR.st.don, DMR.axisSt, color_val_st)
# DMR.st.man
# 
# tiff("./sexDiff/plots/DMR_sexStManPlot.tiff", units="in", width = 7, height = 4, res = 600)
# DMR.st.man
#dev.off()

DMS.st.all.man <- plot.manhat(DMS.st.all.don, DMS.axisSt.all, color_val_st)
DMS.st.all.man

tiff("./sexDiff/plots/DMS_sexStManPlot_all.tiff", units="in", width = 7, height = 4, res = 600)
DMS.st.all.man
dev.off()

DMR.st.all.man <- plot.manhat(DMR.st.all.don, DMR.axisSt.all, color_val_st)
DMR.st.all.man

tiff("./sexDiff/plots/DMR_sexStManPlot_all.tiff", units="in", width = 7, height = 4, res = 600)
DMR.st.all.man
dev.off()

## Analyze annotation stats ## 
GTest(DMS.anno.dev, p = cpg.anno.dev/(cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))
GTest(DMR.anno.dev, p = cpg.anno.dev/(cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

GTest(DMS.anno.st, p = cpg.anno.st/(cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))
GTest(DMR.anno.st, p = cpg.anno.st/(cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

#post hoc tests 
DMS.prom.test.dev <- GTest(c(DMS.anno.dev[1], (DMS.anno.dev[2]+DMS.anno.dev[3]+DMS.anno.dev[4])), 
                  p = c(cpg.anno.dev[1], (cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))/
                  (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

DMS.ex.test.dev <- GTest(c(DMS.anno.dev[2], (DMS.anno.dev[1]+DMS.anno.dev[3]+DMS.anno.dev[4])), 
                    p = c(cpg.anno.dev[2], (cpg.anno.dev[1]+cpg.anno.dev[3]+cpg.anno.dev[4]))/
                      (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

DMS.int.test.dev <- GTest(c(DMS.anno.dev[3], (DMS.anno.dev[1]+DMS.anno.dev[2]+DMS.anno.dev[4])), 
                p = c(cpg.anno.dev[3], (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[4]))/
                  (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

DMS.inter.test.dev <- GTest(c(DMS.anno.dev[4], (DMS.anno.dev[1]+DMS.anno.dev[2]+DMS.anno.dev[3])), 
                      p = c(cpg.anno.dev[4], (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]))/
                        (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))


DMS.prom.test.st <- GTest(c(DMS.anno.st[1], (DMS.anno.st[2]+DMS.anno.st[3]+DMS.anno.st[4])), 
                       p = c(cpg.anno.st[1], (cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))/
                         (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

DMS.ex.test.st <- GTest(c(DMS.anno.st[2], (DMS.anno.st[1]+DMS.anno.st[3]+DMS.anno.st[4])), 
                     p = c(cpg.anno.st[2], (cpg.anno.st[1]+cpg.anno.st[3]+cpg.anno.st[4]))/
                       (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

DMS.int.test.st <- GTest(c(DMS.anno.st[3], (DMS.anno.st[1]+DMS.anno.st[2]+DMS.anno.st[4])), 
                      p = c(cpg.anno.st[3], (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[4]))/
                        (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

DMS.inter.test.st <- GTest(c(DMS.anno.st[4], (DMS.anno.st[1]+DMS.anno.st[2]+DMS.anno.st[3])), 
                        p = c(cpg.anno.st[4], (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]))/
                          (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

#DMRs
DMR.prom.test.dev <- GTest(c(DMR.anno.dev[1], (DMR.anno.dev[2]+DMR.anno.dev[3]+DMR.anno.dev[4])), 
                           p = c(cpg.anno.dev[1], (cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))/
                             (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

DMR.ex.test.dev <- GTest(c(DMR.anno.dev[2], (DMR.anno.dev[1]+DMR.anno.dev[3]+DMR.anno.dev[4])), 
                         p = c(cpg.anno.dev[2], (cpg.anno.dev[1]+cpg.anno.dev[3]+cpg.anno.dev[4]))/
                           (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

DMR.int.test.dev <- GTest(c(DMR.anno.dev[3], (DMR.anno.dev[1]+DMR.anno.dev[2]+DMR.anno.dev[4])), 
                          p = c(cpg.anno.dev[3], (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[4]))/
                            (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))

DMR.inter.test.dev <- GTest(c(DMR.anno.dev[4], (DMR.anno.dev[1]+DMR.anno.dev[2]+DMR.anno.dev[3])), 
                            p = c(cpg.anno.dev[4], (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]))/
                              (cpg.anno.dev[1]+cpg.anno.dev[2]+cpg.anno.dev[3]+cpg.anno.dev[4]))


DMR.prom.test.st <- GTest(c(DMR.anno.st[1], (DMR.anno.st[2]+DMR.anno.st[3]+DMR.anno.st[4])), 
                          p = c(cpg.anno.st[1], (cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))/
                            (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

DMR.ex.test.st <- GTest(c(DMR.anno.st[2], (DMR.anno.st[1]+DMR.anno.st[3]+DMR.anno.st[4])), 
                        p = c(cpg.anno.st[2], (cpg.anno.st[1]+cpg.anno.st[3]+cpg.anno.st[4]))/
                          (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

DMR.int.test.st <- GTest(c(DMR.anno.st[3], (DMR.anno.st[1]+DMR.anno.st[2]+DMR.anno.st[4])), 
                         p = c(cpg.anno.st[3], (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[4]))/
                           (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

DMR.inter.test.st <- GTest(c(DMR.anno.st[4], (DMR.anno.st[1]+DMR.anno.st[2]+DMR.anno.st[3])), 
                           p = c(cpg.anno.st[4], (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]))/
                             (cpg.anno.st[1]+cpg.anno.st[2]+cpg.anno.st[3]+cpg.anno.st[4]))

#adjsut p values 
DMS.dev.gtest.p <- c(DMS.prom.test.dev$p.value,DMS.ex.test.dev$p.value,DMS.int.test.dev$p.value,DMS.inter.test.dev$p.value)
DMS.dev.gtest.p <- p.adjust(DMS.dev.gtest.p, method = "hommel")

DMS.st.gtest.p <- c(DMS.prom.test.st$p.value,DMS.ex.test.st$p.value,DMS.int.test.st$p.value,DMS.inter.test.st$p.value)
DMS.st.gtest.p <- p.adjust(DMS.st.gtest.p, method = "hommel")

DMR.dev.gtest.p <- c(DMR.prom.test.dev$p.value,DMR.ex.test.dev$p.value,DMR.int.test.dev$p.value,DMR.inter.test.dev$p.value)
DMR.dev.gtest.p <- p.adjust(DMR.dev.gtest.p, method = "hommel")

DMR.st.gtest.p <- c(DMR.prom.test.st$p.value,DMR.ex.test.st$p.value,DMR.int.test.st$p.value,DMR.inter.test.st$p.value)
DMR.st.gtest.p <- p.adjust(DMR.st.gtest.p, method = "hommel")

#format data for plotting 
DMS.anno.data.dev <- data.frame(DMS.anno.dev)
DMS.anno.data.dev$dist <- "DMS"
colnames(DMS.anno.data.dev)[1] <- "count"
DMS.anno.data.dev <- rownames_to_column(DMS.anno.data.dev, "feature")

DMR.anno.data.dev <- data.frame(DMR.anno.dev)
DMR.anno.data.dev$dist <- "DMR"
colnames(DMR.anno.data.dev)[1] <- "count"
DMR.anno.data.dev <- rownames_to_column(DMR.anno.data.dev, "feature")

cpg.anno.data.dev <- data.frame(cpg.anno.dev)
cpg.anno.data.dev$dist <- "CpG"
colnames(cpg.anno.data.dev)[1] <- "count"
cpg.anno.data.dev <- rownames_to_column(cpg.anno.data.dev, "feature")

dev.plot.data <- rbind(DMS.anno.data.dev, DMR.anno.data.dev, cpg.anno.data.dev)

DMS.anno.data.st <- data.frame(DMS.anno.st)
DMS.anno.data.st$dist <- "DMS"
colnames(DMS.anno.data.st)[1] <- "count"
DMS.anno.data.st <- rownames_to_column(DMS.anno.data.st, "feature")

DMR.anno.data.st <- data.frame(DMR.anno.st)
DMR.anno.data.st$dist <- "DMR"
colnames(DMR.anno.data.st)[1] <- "count"
DMR.anno.data.st <- rownames_to_column(DMR.anno.data.st, "feature")

cpg.anno.data.st <- data.frame(cpg.anno.st)
cpg.anno.data.st$dist <- "CpG"
colnames(cpg.anno.data.st)[1] <- "count"
cpg.anno.data.st <- rownames_to_column(cpg.anno.data.st, "feature")

st.plot.data <- rbind(DMS.anno.data.st, DMR.anno.data.st, cpg.anno.data.st)

dev.bar.plot.func <- function(data){
          ggplot(data, aes(x = dist, y = count, fill = feature)) +
            geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) + theme_bw() + xlab(NULL) +
            ylab("Percent") + scale_fill_manual(values = c("exon" = "lightgreen", "intron" = "skyblue", 
                                                           "intergenic" = "orange", "promoter" = "yellow"),
                                                name = "Feature") + 
            theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14)) 

}

dev.bar.plot.func <- function(data){
  ggplot(data, aes(x = dist, y = count, fill = feature)) +
    geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) + theme_bw() + xlab(NULL) +
    ylab("Percent") + scale_fill_manual(values = c("exon" = "lightgreen", "intron" = "skyblue", 
                                                   "intergenic" = "orange", "promoter" = "yellow"),
                                        name = "Feature") + 
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) 
  
}

st.bar.plot.func <- function(data){
  ggplot(data, aes(x = dist, y = count, fill = feature)) +
    geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) + theme_bw() + xlab(NULL) +
    ylab("Percent") + scale_fill_manual(values = c("exon" = "#B8DE29FF", "intron" = "#1F968BFF", 
                                                   "intergenic" = "#39568CFF", "promoter" = "#481567FF"),
                                        name = "Feature") + 
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) 
  
}

dev.barplot <- dev.bar.plot.func(dev.plot.data)
dev.barplot

# tiff("./sexDiff/plots/devBarPlot.tiff", units="in", width = 4, height = 5, res = 600)
# dev.barplot
# dev.off()
# 
# st.barplot <- st.bar.plot.func(st.plot.data)
# st.barplot
# 
# tiff("./sexDiff/plots/stBarPlot.tiff", units="in", width = 4, height = 5, res = 600)
# st.barplot
# dev.off()




  