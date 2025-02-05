##################################################################
### Goal: Clean behavioural data and prepare it for analysis 
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
openfield.data <- read_csv("./dev_exp/data/raw/openfield.csv")
shoaling.data <- read_csv("./dev_exp/data/raw/shoaling.csv")
time.data <- read_csv("./dev_exp/data/raw/time.csv")

## Open Field data ##
#merge time data 
openfield.data <- merge(openfield.data, time.data, by="ID")

#export data
write_csv(openfield.data, file("./dev_exp/data/clean/clean_openfield.csv"))

## Shoaling data ##
#change cue labels to match open field data 
shoaling.data$cue <- gsub('alarm', 'ac', shoaling.data$cue)
shoaling.data$cue <- gsub('control', 'c', shoaling.data$cue)

#change missing latency values to 300s (trial length) then convert column to numeric
shoaling.data <- as.data.frame(lapply(shoaling.data, function(y) gsub("-", "300", y)))

#convert to numeric
shoaling.data[,6:29] <- sapply(shoaling.data[,6:29], as.numeric)

#merge time data 
shoaling.data <- merge(shoaling.data, time.data, by="ID")

#export data
write_csv(shoaling.data, file(here("gup_cue_exp", "data", "clean", "clean_shoaling.csv")))

