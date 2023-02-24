#################################
### Goal: Combine all metadata
### Author: Janay Fox
### R script 
#################################

## Set up ##
#install packages
#install.packages("here")
#install.packages("tidyverse")
#install.packages("dplyr")

#load packages
library(here)
library(tidyverse)
library(dplyr)
library(tidyr)

#read in data files
of.data <- read_csv(here("gup_cue_exp", "data", "metadata", "of_metadata.csv"))
sh.data <- read_csv(here("gup_cue_exp", "data", "metadata", "sh_metadata.csv"))

cue.subdiff.5min.data <- read_csv(here("gup_cue_exp", "data", "metadata", "cue_diffsub_fivemin_metadata.csv"))
cue.subdiff.2min.data <- read_csv(here("gup_cue_exp", "data", "metadata", "cue_diffsub_twomin_metadata.csv"))
cue.5min.data <- read_csv(here("gup_cue_exp", "data", "metadata", "cue_fivemin_metadata.csv"))
cue.2min.data <- read_csv(here("gup_cue_exp", "data", "metadata", "cue_twomin_metadata.csv"))
cue.info.data <- read_csv(here("gup_cue_exp", "data", "metadata", "cue_collectioninfo_metadata.csv"))

## Short-term exposure fish ##
#convert data to proportions 
cue.2min.data$sig_prop <- cue.2min.data$sig_s/cue.2min.data$length_obs_s
cue.2min.data$chasing_prop <- cue.2min.data$chasing_s/cue.2min.data$length_obs_s
cue.2min.data$matbehav_prop <- cue.2min.data$matbehav_s/cue.2min.data$length_obs_s
cue.2min.data$bottom_prop <- cue.2min.data$bottom_s/cue.2min.data$length_obs_s

cue.5min.data$sig_prop <- cue.5min.data$sig_s/cue.5min.data$length_obs_s
cue.5min.data$chasing_prop <- cue.5min.data$chasing_s/cue.5min.data$length_obs_s
cue.5min.data$matbehav_prop <- cue.5min.data$matbehav_s/cue.5min.data$length_obs_s
cue.5min.data$bottom_prop <- cue.5min.data$bottom_s/cue.5min.data$length_obs_s

#remove unneccassary columns 
cue.2min.data <- cue.2min.data[-c(6:17, 19:20,23)]
cue.5min.data <- cue.5min.data[-c(6:17, 19:20,23)]

#convert to wide format
cue.2min.data <- cue.2min.data %>% pivot_wider(names_from = time, values_from = 
                              c(dash_nb, subuse_prop, frozen_prop, sig_prop, 
                              chasing_prop, matbehav_prop, bottom_prop))

cue.5min.data <- cue.5min.data %>% pivot_wider(names_from = time, values_from = 
                                                 c(dash_nb, subuse_prop, frozen_prop, sig_prop, 
                                                   chasing_prop, matbehav_prop, bottom_prop))

#rename columns to specify trial length 
cue.2min.data <- cue.2min.data %>% rename("dash_nb_2min_before" = "dash_nb_before",
                                          "dash_nb_2min_after" ="dash_nb_after",
                                          "subuse_prop_2min_before" = "subuse_prop_before",
                                          "subuse_prop_2min_after" = "subuse_prop_after" ,
                                          "frozen_prop_2min_before" = "frozen_prop_before",
                                          "frozen_prop_2min_after" = "frozen_prop_after",
                                          "sig_prop_2min_before" = "sig_prop_before",
                                          "sig_prop_2min_after" = "sig_prop_after",
                                          "chasing_prop_2min_before" = "chasing_prop_before",
                                          "chasing_prop_2min_after" = "chasing_prop_after", 
                                          "matbehav_prop_2min_before" ="matbehav_prop_before",
                                          "matbehav_prop_2min_after" = "matbehav_prop_after",
                                          "bottom_prop_2min_before" = "bottom_prop_before",
                                          "bottom_prop_2min_after" = "bottom_prop_after")

cue.5min.data <- cue.5min.data %>% rename("dash_nb_5min_before" = "dash_nb_before",
                                          "dash_nb_5min_after" ="dash_nb_after",
                                          "subuse_prop_5min_before" = "subuse_prop_before",
                                          "subuse_prop_5min_after" = "subuse_prop_after" ,
                                          "frozen_prop_5min_before" = "frozen_prop_before",
                                          "frozen_prop_5min_after" = "frozen_prop_after",
                                          "sig_prop_5min_before" = "sig_prop_before",
                                          "sig_prop_5min_after" = "sig_prop_after",
                                          "chasing_prop_5min_before" = "chasing_prop_before",
                                          "chasing_prop_5min_after" = "chasing_prop_after", 
                                          "matbehav_prop_5min_before" ="matbehav_prop_before",
                                          "matbehav_prop_5min_after" = "matbehav_prop_after",
                                          "bottom_prop_5min_before" = "bottom_prop_before",
                                          "bottom_prop_5min_after" = "bottom_prop_after")

cue.subdiff.2min.data <- cue.subdiff.2min.data %>% rename("diff_sub_use_2min" = "diff_sub_use")
cue.subdiff.5min.data <- cue.subdiff.5min.data %>% rename("diff_sub_use_5min" = "diff_sub_use")

#merge datasets
cue.metadata <- merge(cue.2min.data, cue.5min.data, by=c("fish_ID", "sex", "cue", "tank"))
cue.metadata <- merge(cue.metadata,cue.info.data, by=c("fish_ID"))
cue.metadata <- merge(cue.metadata, cue.subdiff.2min.data, by=c("fish_ID", "sex", "cue", "tank"))
cue.metadata <- merge(cue.metadata, cue.subdiff.5min.data, by=c("fish_ID", "sex", "cue", "tank"))

#export data
write_csv(cue.metadata, file(here("gup_cue_exp", "data", "clean", "clean_shortterm_metadata.csv")))

## Developmental Fish ##
#merge dataframes
dev.metadata <- merge(of.data, sh.data, by=c("ID", "sex", "cue", "tank"))

#edit column name to be consistent with other dataset 
dev.metadata <- dev.metadata %>% rename("fish_ID" = "ID")

#edit column values to be consistent with other dataset 
dev.metadata$sex[dev.metadata$sex== "f"] <- "female"
dev.metadata$sex[dev.metadata$sex== "m"] <- "male"

dev.metadata$cue[dev.metadata$cue== "ac"] <- "alarm"
dev.metadata$cue[dev.metadata$cue== "c"] <- "control"

dev.metadata$fish_ID <- paste("D", dev.metadata$fish_ID, sep = "")
dev.metadata$tank <- paste("D", dev.metadata$tank, sep = "")

#remove unnecassry columns 
dev.metadata <- dev.metadata[-c(6:7,9:10, 13,16,23:35)]

## Export data ##
write_csv(dev.metadata, file(here("gup_cue_exp", "data", "clean", "clean_developmental_metadata.csv")))
