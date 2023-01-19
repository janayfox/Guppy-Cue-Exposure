############################################
### Goal: Plot and analyze weight/length data
### Author: Janay Fox
### R script
############################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("here")
#install.packages("forcats")
#install.packages("reshape2")
#install.packages("factoextra")
#install.packages("tidyverse")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("ggpubr")
#install.packages("MuMIn")
#install.packages("DHARMa")

#load packaes
library(ggplot2)
library(dplyr)
library(here)
library(forcats)
library(reshape2)
library(factoextra)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggpubr)
library(MuMIn)
library(DHARMa)
library(ggstatsplot)

#read data files
data <- read_csv(here("gup_cue_exp", "data", "clean", "dev_size.csv")) #import weight data 

## Check for impact of cue on weight, length, condition 
#calculate condition 
data$condition <- (100 * data$weight_g) / (data$tl_cm)^3

#plot data
length.box <- data %>% ggbetweenstats(x = cue, y = sl_cm,  
                                          messages = FALSE, results.subtitle = FALSE, pairwise.comparisons = TRUE,
                                          p.adjust.method = "bonferroni", centrality.plotting = FALSE) + 
  xlab("Cue") + ylab("Standard Length (cm)") + theme_bw() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 15), legend.position = "none") + 
  scale_x_discrete(labels = c("ac" = "Alarm Cue", "c" = "Control")) 

weight.box <- data %>% ggbetweenstats(x = cue, y = weight_g,  
                                      messages = FALSE, results.subtitle = FALSE, pairwise.comparisons = TRUE,
                                      p.adjust.method = "bonferroni", centrality.plotting = FALSE) + 
  xlab("Cue") + ylab("Weight (g)") + theme_bw() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 15), legend.position = "none") + 
  scale_x_discrete(labels = c("ac" = "Alarm Cue", "c" = "Control")) 

cond.box <- data %>% ggbetweenstats(x = cue, y = condition,  
                                      messages = FALSE, results.subtitle = FALSE, pairwise.comparisons = TRUE,
                                      p.adjust.method = "bonferroni", centrality.plotting = FALSE) + 
  xlab("Cue") + ylab("Condition") + theme_bw() +
  theme(axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 15), legend.position = "none") + 
  scale_x_discrete(labels = c("ac" = "Alarm Cue", "c" = "Control")) 

#no significant differences bewteen cues
