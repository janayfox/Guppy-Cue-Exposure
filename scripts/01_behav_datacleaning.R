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
library(here)
library(tidyverse)
library(dplyr)

#read in data files
openfield.data <- read_csv(here("gup_cue_exp", "data", "raw", "openfield.csv"))
feeding.data <- read_csv(here("gup_cue_exp","data", "raw", "feeding.csv"))
shoaling.data <- read_csv(here("gup_cue_exp","data", "raw", "shoaling.csv"))
SL.data <- read_csv(here("gup_cue_exp","data", "raw", "sociallearning.csv"))
SL.lat.data <- read_csv(here("gup_cue_exp","data", "raw", "SL_lat.csv"))
SL.dem.data <- read_csv(here("gup_cue_exp","data", "raw", "sl_dem.csv"))
time.data <- read_csv(here("gup_cue_exp", "data", "raw", "time.csv"))

## Open Field data ##
#merge time data 
openfield.data <- merge(openfield.data, time.data, by="ID")

#export data
write_csv(openfield.data, file(here("gup_cue_exp", "data", "clean", "clean_openfield.csv")))

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

## Feeding data ##
#add on sex, tank and cue columns 
feeding.data <- merge(feeding.data, openfield.data[1:4], by="ID")

#merge time data 
feeding.data <- merge(feeding.data, time.data, by="ID")

#export data
write_csv(feeding.data, file(here("gup_cue_exp", "data", "clean", "clean_feeding.csv")))

### Social Learning data ###
#combine bins to get full trial values
SL_full.data <- SL.data %>% group_by(fish_ID) %>% summarize(in_L_s = sum(in_L_s), 
                                                            in_L_nb = sum(in_L_nb),
                                                            in_R_s = sum(in_R_s),
                                                            in_R_nb = sum(in_R_nb),
                                                            nr_L_s = sum(nr_L_s),
                                                            nr_L_nb = sum(nr_L_nb),
                                                            nr_R_s = sum(nr_R_s),
                                                            nr_R_nb = sum(nr_R_nb),
                                                            far_L_s = sum(far_L_s),
                                                            far_L_nb = sum(far_L_nb),
                                                            far_R_s = sum(far_R_s),
                                                            far_R_nb = sum(far_R_nb),
                                                            peck_L_nb = sum(peck_L_nb),
                                                            peck_R_nb = sum(peck_R_nb))
#subset bin 1 for 5 min trials
SL_5min.data <- subset(SL.data, Bin == 1)
SL_5min.data <- SL_5min.data[c(-1,-2)] #drop obs id and bin columns

#calculate latency 
SL.lat.data$lat_in_L <- (SL.lat.data$time_in_L - SL.lat.data$video_start)
SL.lat.data$lat_in_R <- (SL.lat.data$time_in_R - SL.lat.data$video_start)
SL.lat.data$lat_nr_L <- (SL.lat.data$time_near_L - SL.lat.data$video_start)
SL.lat.data$lat_nr_R <- (SL.lat.data$time_near_R - SL.lat.data$video_start)
SL.lat.data$lat_far_L <- (SL.lat.data$time_far_L - SL.lat.data$video_start)
SL.lat.data$lat_far_R <- (SL.lat.data$time_far_R - SL.lat.data$video_start)

#combine latency values onto data frame 
SL_full.data <- merge(SL_full.data, SL.lat.data[,c(1,9:14)], by = "fish_ID")

SL_5min.data <- merge(SL_5min.data, SL.lat.data[,c(1,9:14)], by = "fish_ID")

#change ID column name to match other dataframes 
colnames(SL_full.data)[1] <- "ID"
colnames(SL_5min.data)[1] <- "ID"

#add on sex, tank and cue columns 
SL_full.data <- merge(SL_full.data, openfield.data[1:4], by ="ID")
SL_5min.data <- merge(SL_5min.data, openfield.data[1:4], by ="ID")

#add on demonstrated feeder data
SL_full.data <- merge(SL_full.data, SL.dem.data, by ="ID")
SL_5min.data <- merge(SL_5min.data, SL.dem.data, by ="ID")

#make columns for time in dem/nondem feeders
SL_full.data$in_dem_s <- NA
SL_full.data$in_dem_lat_s <- NA
SL_full.data$in_nondem_s <- NA
SL_full.data$in_nondem_lat_s <- NA
SL_full.data$nr_dem_s <- NA
SL_full.data$nr_dem_lat_s <- NA
SL_full.data$nr_nondem_s <- NA
SL_full.data$nr_nondem_lat_s <- NA
SL_full.data$far_dem_s <- NA
SL_full.data$far_dem_lat_s <- NA
SL_full.data$far_nondem_s <- NA
SL_full.data$far_nondem_lat_s <- NA

SL_5min.data$in_dem_s <- NA
SL_5min.data$in_dem_lat_s <- NA
SL_5min.data$in_nondem_s <- NA
SL_5min.data$in_nondem_lat_s <- NA
SL_5min.data$nr_dem_s <- NA
SL_5min.data$nr_dem_lat_s <- NA
SL_5min.data$nr_nondem_s <- NA
SL_5min.data$nr_nondem_lat_s <- NA
SL_5min.data$far_dem_s <- NA
SL_5min.data$far_dem_lat_s <- NA
SL_5min.data$far_nondem_s <- NA
SL_5min.data$far_nondem_lat_s <- NA

#create function that takes x,y,z and converts columns of L/R feeder to dem/nondem feeding
# x = the column you want to fill in 
# y = the R/hor column for dem columns and L/ver for nondem columns
# z = the L/ver column for dem columns and R/ver for nondem columns
conv_col_full <- function(x,y,z) { for (i in 1:79){
  if (SL_full.data[i,25] == "hor"){
  SL_full.data[i,x] = SL_full.data[i,y]}
    else if (SL_full.data[i,25] == "ver"){
      SL_full.data[i,x] = SL_full.data[i,z]}
}
  return(SL_full.data)}

conv_col_5min <- function(x,y,z) { for (i in 1:79){
  if (SL_5min.data[i,25] == "hor"){
    SL_5min.data[i,x] = SL_5min.data[i,y]}
  else if (SL_5min.data[i,25] == "ver"){
    SL_5min.data[i,x] = SL_5min.data[i,z]}
}
  return(SL_5min.data)}

#populate columns with proper values depending on demonstrated feeder 
SL_full.data <- conv_col_full(26,4,2)
SL_full.data <- conv_col_full(27,17,16)
SL_full.data <- conv_col_full(28,2,4)
SL_full.data <- conv_col_full(29,16,17)
SL_full.data <- conv_col_full(30,8,6)
SL_full.data <- conv_col_full(31,19,18)
SL_full.data <- conv_col_full(32,6,8)
SL_full.data <- conv_col_full(33,18,19)
SL_full.data <- conv_col_full(34,12,10)
SL_full.data <- conv_col_full(35,21,20)
SL_full.data <- conv_col_full(36,10,12)
SL_full.data <- conv_col_full(37,20,21)

SL_5min.data <- conv_col_5min(26,4,2)
SL_5min.data <- conv_col_5min(27,17,16)
SL_5min.data <- conv_col_5min(28,2,4)
SL_5min.data <- conv_col_5min(29,16,17)
SL_5min.data <- conv_col_5min(30,8,6)
SL_5min.data <- conv_col_5min(31,19,18)
SL_5min.data <- conv_col_5min(32,6,8)
SL_5min.data <- conv_col_5min(33,18,19)
SL_5min.data <- conv_col_5min(34,12,10)
SL_5min.data <- conv_col_5min(35,21,20)
SL_5min.data <- conv_col_5min(36,10,12)
SL_5min.data <- conv_col_5min(37,20,21)

#change NAs to full time of trial 
SL_full.data[is.na(SL_full.data)] = 600
SL_5min.data[is.na(SL_5min.data)] = 300

#replace latency values over 300 with 300 in 5min trial 
SL_5min.data$in_dem_lat_s <- replace(SL_5min.data$in_dem_lat_s, SL_5min.data$in_dem_lat_s>300, 300)
SL_5min.data$in_nondem_lat_s <- replace(SL_5min.data$in_nondem_lat_s, SL_5min.data$in_nondem_lat_s>300, 300)
SL_5min.data$nr_dem_lat_s <- replace(SL_5min.data$nr_dem_lat_s, SL_5min.data$nr_dem_lat_s>300, 300)
SL_5min.data$nr_nondem_lat_s <- replace(SL_5min.data$nr_nondem_lat_s, SL_5min.data$nr_nondem_lat_s>300, 300)
SL_5min.data$far_dem_lat_s <- replace(SL_5min.data$far_dem_lat_s, SL_5min.data$far_dem_lat_s>300, 300)
SL_5min.data$far_nondem_lat_s <- replace(SL_5min.data$far_nondem_lat_s, SL_5min.data$far_nondem_lat_s>300, 300)

#merge time data 
SL_full.data <- merge(SL_full.data, time.data, by="ID")
SL_5min.data <- merge(SL_5min.data, time.data, by="ID")

#export data
write_csv(SL_full.data, file(here("gup_cue_exp", "data", "clean", "clean_full_sl.csv")))
write_csv(SL_5min.data, file(here("gup_cue_exp", "data", "clean", "clean_5min_sl.csv")))
