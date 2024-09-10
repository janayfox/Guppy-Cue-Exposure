##############################################################
### Goal: Get stats for alignment and number of reads
### Author: Janay Fox
### R script
#############################################################

## Set up ##
#install packages
#install.packages("here")
#install.packages("tidyverse")
#install.packages("dplyr")

#load packages
library(here)
library(tidyverse)
library(dplyr)

#read in data files
align.dev.data <- read_csv(here("gup_cue_exp", "data", "alignment_res", "alignment_meth_dev.csv"))
align.st.data <- read_csv(here("gup_cue_exp", "data", "alignment_res", "alignment_meth_shortterm.csv"))

## Stats ##
#calculate means and standard devs
#reads
mean(align.dev.data$Reads)
sd(align.dev.data$Reads)

mean(align.st.data$Reads)
sd(align.st.data$Reads)

#mapping efficiency 
mean(align.dev.data$Perc_Mapping_Efficiency)
sd(align.dev.data$Perc_Mapping_Efficiency)

mean(align.st.data$Perc_Mapping_Efficiency)
sd(align.st.data$Perc_Mapping_Efficiency)

#Number of CpGs
mean(align.dev.data$No_CpGs)
sd(align.dev.data$No_CpGs)

mean(align.st.data$No_CpGs)
sd(align.st.data$No_CpGs)

#percent methylated CpGs 
mean(align.dev.data$Perc_Meth_CpGs)
sd(align.dev.data$Perc_Meth_CpGs)

#in each group
tapply(align.dev.data$Perc_Meth_CpGs, align.dev.data$Treatment, mean)
tapply(align.dev.data$Perc_Meth_CpGs, align.dev.data$Treatment, sd)

mean(align.st.data$Perc_Meth_CpGs)
sd(align.st.data$Perc_Meth_CpGs)

#percent methylated CHGs 
mean(align.dev.data$Perc_Meth_CHGs)
sd(align.dev.data$Perc_Meth_CHGs)

mean(align.st.data$Perc_Meth_CHGs)
sd(align.st.data$Perc_Meth_CHGs)

#percent methylated CHGs 
mean(align.dev.data$Perc_Meth_CHHs)
sd(align.dev.data$Perc_Meth_CHHs)

mean(align.st.data$Perc_Meth_CHHs)
sd(align.st.data$Perc_Meth_CHHs)

