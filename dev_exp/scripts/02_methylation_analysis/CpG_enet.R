#library(ggplot2)
library(dplyr)
library(methylKit)
#library(lme4)
#library(car)
#library(plsRglm)
#library(mixOmics)
library(glmnet)
library(caret)

#read behaviour data 
of.data <- read.csv("./clean_openfield.csv")
sh.data <- read.csv("./clean_shoaling.csv")

#change ID name to match 
of.data$ID <- paste0("D", of.data$ID)
sh.data$ID <- paste0("D", sh.data$ID)

#create new boldness variable 
of.data$boldness_score <- of.data$time_shelter_s + of.data$frozen_s

#create new exploration variable
cumDur_data <- of.data[,73:100]

unique_data_cumDur <- cumDur_data %>%
  rowwise() %>%
  mutate(unique_squares = sum(c_across(starts_with("S")) > 3)) %>%
  ungroup()

of.data$uni_sq <- unique_data_cumDur$unique_squares

#create new shoaling variable
sh.data$total_shoal_s <- sh.data$tight_shoal_s + sh.data$lse_shoal_s
sh.data$total_emp_s <- sh.data$tight_emp_s + sh.data$lse_emp_s
sh.data$total_diff_s <- sh.data$total_shoal_s - sh.data$total_emp_s

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

#make function to transpose data and remove NAs
prep_data <- function(df){
  #remove NAs
  df <- na.omit(df)
  
  #transpose data
  t.df <- t(df)
  
  
  return(as.data.frame(t.df))
}

#function to get best results from caret fit
get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

#Prepare data ##
#get data for all CpGs
CpG_meth_fem <- readRDS("./DMSmeth_fem_5X.RDS")
#CpG_meth_mal <- readRDS("./DMSmeth_mal_5X.RDS")

#calculate percent methylation
CpG_fem_percMeth <- as.data.frame(percMethylation(CpG_meth_fem))
#CpG_mal_percMeth <- as.data.frame(percMethylation(CpG_meth_mal))

#add site names
CpG_meth_fem <- nameSite(CpG_meth_fem)
#CpG_meth_mal <- nameSite(CpG_meth_mal)

#add rownames to perc meth matrix
rownames(CpG_fem_percMeth) <- CpG_meth_fem$site_name
#rownames(CpG_mal_percMeth) <- CpG_meth_mal$site_name

#remove extra file
rm(CpG_meth_fem)
#rm(CpG_meth_mal)

#transpose
CpG_fem_percMeth <- prep_data(CpG_fem_percMeth)
#CpG_mal_percMeth <- prep_data(CpG_mal_percMeth)

#make ID column
CpG_fem_percMeth$ID <- rownames(CpG_fem_percMeth)
#CpG_mal_percMeth$ID <- rownames(CpG_mal_percMeth)

#create behavioural dataframes
bold.data <- of.data[,c("ID","boldness_score")]
exp.data <- of.data[,c("ID", "uni_sq")]
act.data <- of.data[,c("ID", "dist_cm")]
shoal.data <- sh.data[,c("ID", "total_diff_s")]

## Run Elastic Net ## 
#set parameters
set.seed(69)
cv.10 = trainControl(method = "cv", number = 10, verboseIter = TRUE, savePredictions = "final")

#function to run enet without splitting data
enet.caret.function.v2 <- function(behav.data, meth.data, behav.col, fam){
  #merge data then remove ID
  merged.data <- merge(behav.data, meth.data, by = "ID")
  merged.data <- merged.data[, !(names(merged.data) %in% "ID")]
  
  #train
  enet.model <- train(
    as.formula(paste(behav.col, "~ .")),
    data = merged.data,
    method = "glmnet",
    trControl = cv.10,
    tuneLength = 10,
    family = fam,
    preProcess = c("center", "scale")
  )
  
  #print best model
  print(get_best_result(enet.model))
  
  #get variable importance 
  var.importance <- varImp(enet.model, scale = FALSE)
  
  #return 
  return(list(model = enet.model, variable_importance = var.importance))
}

#CpG
fem.bold.DMS <- enet.caret.function(bold.data, CpG_fem_percMeth, "boldness_score", "gaussian")
#mal.bold.DMS <- enet.caret.function(bold.data, CpG_mal_percMeth, "boldness_score", "gaussian")

fem.exp.DMS <- enet.caret.function(exp.data, CpG_fem_percMeth, "uni_sq", "poisson")
#mal.exp.DMS <- enet.caret.function(exp.data, CpG_mal_percMeth, "uni_sq", "poisson")

fem.act.DMS <- enet.caret.function(act.data, CpG_fem_percMeth, "dist_cm", "gaussian")
#mal.act.DMS <- enet.caret.function(act.data, CpG_mal_percMeth, "dist_cm", "gaussian")

fem.shoal.DMS <- enet.caret.function(shoal.data, CpG_fem_percMeth, "total_diff_s", "gaussian")
#mal.shoal.DMS <- enet.caret.function(shoal.data, CpG_mal_percMeth, "total_diff_s","gaussian")

#save RDS
saveRDS(fem.bold.DMS, "./fem_bold_CpG_enet.RDS")
#saveRDS(mal.bold.DMS, "./mal_bold_CpG_enet.RDS")

saveRDS(fem.exp.DMS, "./fem_exp_CpG_enet.RDS")
#saveRDS(mal.exp.DMS, "./mal_exp_CpG_enet.RDS")

saveRDS(fem.act.DMS, "./fem_act_CpG_enet.RDS")
#saveRDS(mal.act.DMS, "./mal_act_CpG_enet.RDS")

saveRDS(fem.shoal.DMS, "./fem_shoal_CpG_enet.RDS")
#saveRDS(mal.shoal.DMS, "./mal_shoal_CpG_enet.RDS")

