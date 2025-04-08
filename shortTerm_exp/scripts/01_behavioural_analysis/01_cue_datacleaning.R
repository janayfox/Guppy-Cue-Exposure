##################################################################
### Goal: Clean cue response data and prepare it for analysis 
### Author: Janay Fox
### R script
##################################################################

## Set up ##
#install packages
#install.packages("here")
#install.packages("tidyverse")
#install.packages("dplyr")

#load packages
library(tidyverse)
library(dplyr)

#read in data files
cue.data <- read_csv("./shortTerm_exp/data/raw/cueresponse.csv")
cue.tank.data <- read_csv("./shortTerm_exp/data/raw/cuetanks.csv")

## Data cleanup ##
#create new column for video bins based off observation ID (contains bin number in last 3 characters)
cue.data$bin <- str_sub(cue.data$obs_id, -3)
cue.data$bin <- as.numeric(cue.data$bin)

#trials that had video broken up into multiple parts need bin number manually adjusted
cue.data$bin[cue.data$obs_id == "C02V3000"] <- 2
cue.data$bin[cue.data$obs_id == "C02V3001"] <- 3
cue.data$bin[cue.data$obs_id == "C02V3002"] <- 4
cue.data$bin[cue.data$obs_id == "C02V3003"] <- 5
cue.data$bin[cue.data$obs_id == "C02V3004"] <- 6
cue.data$bin[cue.data$obs_id == "C02V3005"] <- 7
cue.data$bin[cue.data$obs_id == "C02V3006"] <- 8
cue.data$bin[cue.data$obs_id == "C02V3007"] <- 9

cue.data$bin[cue.data$obs_id == "C04V3000"] <- 8
cue.data$bin[cue.data$obs_id == "C04V3001"] <- 9
cue.data$bin[cue.data$obs_id == "C04V3002"] <- 10

#remove bin from end of observation ID
cue.data$obs_id <- str_sub(cue.data$obs_id,1, nchar(cue.data$obs_id)-3)

#convert character columns to factors
cue.data[sapply(cue.data, is.character)] <- lapply(cue.data[sapply(cue.data, is.character)], as.factor)

#merge with tank dataframe 
cue.data <- merge(cue.data, cue.tank.data, by=c("obs_id","sex"))

#subset into before and after to do different operations on 
before.cue.data <- cue.data %>% subset(time == "before")
after.cue.data <- cue.data %>% subset(time == "after")

## Combine bins for different time frames - 2 min and 5 min ##
#before cue
#subset into last 5 min then aggregate by fish ID and time
fivemin.before.cue.data <- subset(before.cue.data, bin >= 5 & bin <= 9) %>% group_by(fish_ID, sex, time, cue, tank) %>% summarize(length_obs_s = sum(length_obs_s),
                                                                                                               sig_s = sum(sig_s), sig_nb = sum(sig_nb), 
                                                                                                               chasing_s = sum(chasing_s),
                                                                                                            feeding_s = sum(feeding_s), frozen_s = sum(frozen_s), 
                                                                                                            Q1_s = sum(Q1_s), Q2_s = sum(Q2_s), Q3_s = sum(Q3_s),
                                                                                                            Q4_s = sum(Q4_s), Q5_s = sum(Q5_s), Q6_s = sum(Q6_s),
                                                                                                           dash_nb = sum(dash_nb))

#after cue
fivemin.after.cue.data <- subset(after.cue.data, bin >= 0 & bin <= 4) %>% group_by(fish_ID, sex, time, cue, tank) %>% summarize(length_obs_s = sum(length_obs_s),
                                                                                                                             sig_s = sum(sig_s), sig_nb = sum(sig_nb), 
                                                                                                                             chasing_s = sum(chasing_s),
                                                                                                                             feeding_s = sum(feeding_s), frozen_s = sum(frozen_s), 
                                                                                                                             Q1_s = sum(Q1_s), Q2_s = sum(Q2_s), Q3_s = sum(Q3_s),
                                                                                                                             Q4_s = sum(Q4_s), Q5_s = sum(Q5_s), Q6_s = sum(Q6_s),
                                                                                                                             dash_nb = sum(dash_nb))

#merge before and after subsets 
fivemin.cue.data <- rbind(fivemin.before.cue.data,fivemin.after.cue.data)

#export data 
write_csv(fivemin.cue.data, "./shortTerm_exp/data/clean/clean_cue_fivemin.csv")
write_csv(cue.data, "./shortTerm_exp/data/clean/all_cue_data.csv")
