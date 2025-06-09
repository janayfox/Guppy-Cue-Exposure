#######################################################
### Goal: Try PLSR with plsRglm
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("methylKit")
#install.packages(lme4)
#install.packages(lmerTest)
#install.packages(car)
#install.packages(plsRglm)
#install.packages(mixOmics)
#install.packages(glmnet)
#install.packages(caret)
#install.packages(r2glmm)
#install.packages(DHARMa)
#install.packages(S4Vectors)
#install.packages(IRanges)
#install.packages(GenomicRanges)
#install.packages(genomation)
#install.packages(selectiveInference)
#install.packages(ggeffects)
#install.packages(ggpubr)

library(ggplot2)
library(dplyr)
library(methylKit)
library(lme4)
library(lmerTest)
library(car)
library(plsRglm)
library(mixOmics)
library(glmnet)
library(caret)
library(r2glmm)
library(DHARMa)
library(S4Vectors)
library(IRanges)
library(GenomicRanges)
library(genomation)
library(selectiveInference)
library(ggeffects)
library(ggpubr)

#read behaviour data 
sh.data <- read.csv("./dev_exp/data/clean/clean_shoaling.csv")

#change ID name to match 
sh.data$ID <- paste0("D", sh.data$ID)

#create new shoaling variable
sh.data$total_shoal_s <- sh.data$tight_shoal_s + sh.data$lse_shoal_s
sh.data$total_emp_s <- sh.data$tight_emp_s + sh.data$lse_emp_s
sh.data$total_diff_s <- sh.data$total_shoal_s - sh.data$total_emp_s

#get just tank and cue data 
tank.data <- sh.data[,c("ID", "tank", "cue")]

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
  
  #return 
  return(enet.model)
}

enet.caret.function.v3 <- function(behav.data, meth.data, behav.col, fam) {
  # Merge data then remove ID
  merged.data <- merge(behav.data, meth.data, by = "ID")
  merged.data <- merged.data[, !(names(merged.data) %in% "ID")]
  
  # Define leave-one-out cross-validation
  loo <- trainControl(method = "LOOCV")
  
  # Train the elastic net model using leave-one-out cross-validation
  enet.model <- train(
    as.formula(paste(behav.col, "~ .")),
    data = merged.data,
    method = "glmnet",
    trControl = loo,       # Changed from cv.10 to loo
    tuneLength = 10,
    family = fam,
    preProcess = c("center", "scale")
  )
  
  # Print the best model
  print(get_best_result(enet.model))
  
  # Return the trained model
  return(enet.model)
}

#make function to get top coefficients 
get.top.coefficients <- function(enet.model, top_n = 5){
  # Get the best lambda value
  best_lambda <- enet.model$bestTune$lambda
  
  # Extract the coefficients at the best lambda from the final model
  best_model <- enet.model$finalModel
  coefficients <- coef(best_model, s = best_model$lambdaOpt)
  
  # Convert sparse matrix to a readable format
  coefficients <- as.matrix(coefficients)
  
  # Remove intercept (assuming it's the first row)
  coefficients <- coefficients[-1, , drop = FALSE]
  
  # Select non-zero coefficients
  non_zero_coefficients <- coefficients[coefficients != 0, , drop = FALSE]
  
  # Create a data frame with variable names and their absolute coefficients
  coef_df <- data.frame(
    variable = rownames(non_zero_coefficients),
    coefficient = non_zero_coefficients[, 1]
  )
  
  # Order by the absolute value of the coefficients
  coef_df <- coef_df[order(-abs(coef_df$coefficient)), ]
  
  # Select the top N coefficients
  top_coefficients <- head(coef_df, top_n)
  
  return(top_coefficients)
}

#function for running linear models
run.lm.func <- function(meth.data,behav.data,sites,behav.col){
  #merge behaviour data with tank data 
  behav.data <- merge(behav.data, tank.data, by = "ID")
  #select only top sites
  data <- meth.data[,colnames(meth.data) %in% sites]
  
  #scale methylation data
  data <- scale(data)
  
  #set rownames 
  rownames(behav.data) <- behav.data$ID
  behav.data$ID <- NULL
  
  #merge with behavioural data
  data <- merge(data,behav.data, by = "row.names")
  data$Row.names <- NULL
  
  # Dynamically create the formula
  predictors <- setdiff(colnames(data), c("Row.names", "tank", behav.col))
  formula_str <- paste(behav.col, "~", paste(predictors, collapse = " + "), "+ (1|tank)")
  formula <- as.formula(formula_str)
  
  print(formula)
  
  #run formula
  model <- lmer(formula, data = data, REML = TRUE)
  
  print(summary(model))
  
  #calculate stats with car package 
  print(car::Anova(model, type = 2))
  
  #calculate partial r2
  print(r2beta(model, partial = TRUE,  method = 'nsj'))
  
  return(model)
}

#function for DHARMA model validation 
val.func <- function(model){
  simOutput_dist <- simulateResiduals(fittedModel = model)
  return(plot(simOutput_dist))
}

#create function that renames chromosomes
renameChr <- function(obj) {
  obj.gr <- as(obj, "GRanges")
  gr.obj.rename <- renameSeqlevels(obj.gr, c(NC_024331.1="LG1", NC_024332.1="LG2",
                                             NC_024333.1="LG3", NC_024334.1="LG4",
                                             NC_024335.1="LG5",NC_024336.1="LG6",
                                             NC_024337.1="LG7", NC_024338.1="LG8",
                                             NC_024339.1="LG9", NC_024340.1="LG10",
                                             NC_024341.1="LG11", NC_024342.1="LG12",
                                             NC_024343.1="LG13", 
                                             NC_024344.1="LG14", NC_024345.1="LG15",
                                             NC_024346.1="LG16", NC_024347.1="LG17",
                                             NC_024348.1="LG18", NC_024349.1="LG19",
                                             NC_024350.1="LG20", NC_024351.1="LG21",
                                             NC_024352.1="LG22", NC_024353.1="LG23"))
  return(gr.obj.rename)
}

#Prepare data ##
#get percent methylation 
DMS_fem_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMS_res/DMSdiffmeth_fem_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMS_res/DMSmeth_fem_5X.RDS")
DMS_fem_percMeth <- prep_data(DMS_fem_percMeth)

DMS_mal_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMS_res/DMSdiffmeth_mal_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMS_res/DMSmeth_mal_5X.RDS")
DMS_mal_percMeth <- prep_data(DMS_mal_percMeth)

#get perc meth for DMRs
DMR_fem_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMR_res/DMRdiffmeth_fem_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMR_res/DMR_tile_meth_fem_5X.RDS")
DMR_fem_percMeth <- prep_data(DMR_fem_percMeth)

DMR_mal_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMR_res/DMRdiffmeth_mal_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMR_res/DMR_tile_meth_mal_5X.RDS")
DMR_mal_percMeth <- prep_data(DMR_mal_percMeth)

#add ID columns
DMS_fem_percMeth$ID <- rownames(DMS_fem_percMeth)
DMS_mal_percMeth$ID <- rownames(DMS_mal_percMeth)

DMR_fem_percMeth$ID <- rownames(DMR_fem_percMeth)
DMR_mal_percMeth$ID <- rownames(DMR_mal_percMeth)

#create behavioural dataframe
shoal.data <- sh.data[,c("ID", "total_diff_s")]

## Run Elastic Net ## 
#set parameters
set.seed(69)
cv.10 = trainControl(method = "repeatedcv", number = 10, repeats = 10, savePredictions = "all", selectionFunction = "oneSE")

#DMS
fem.shoal.DMS <- enet.caret.function.v2(shoal.data, DMS_fem_percMeth, "total_diff_s", "gaussian")
mal.shoal.DMS <- enet.caret.function.v2(shoal.data, DMS_mal_percMeth, "total_diff_s","gaussian")

#DMR
fem.shoal.DMR <- enet.caret.function.v2(shoal.data, DMR_fem_percMeth, "total_diff_s", "gaussian")
mal.shoal.DMR <- enet.caret.function.v2(shoal.data, DMR_mal_percMeth, "total_diff_s","gaussian")

#try running separate for AC and C 
#doing for shoaling model only as this is the only behaviour I retain in the paper
#make lists of AC ind and C ind 
ac.id <- subset(sh.data, cue == "ac")$ID
c.id <- subset(sh.data, cue == "c")$ID

# make ac and c only datasets 
shoal.data.ac <- shoal.data[shoal.data$ID %in% ac.id,]
shoal.data.c <- shoal.data[shoal.data$ID %in% c.id,]
DMS_fem_percMeth.ac <- DMS_fem_percMeth[DMS_fem_percMeth$ID %in% ac.id,]
DMS_fem_percMeth.c <- DMS_fem_percMeth[DMS_fem_percMeth$ID %in% c.id,]
DMS_mal_percMeth.ac <- DMS_mal_percMeth[DMS_mal_percMeth$ID %in% ac.id,]
DMS_mal_percMeth.c <- DMS_mal_percMeth[DMS_mal_percMeth$ID %in% c.id,]

DMR_fem_percMeth.ac <- DMR_fem_percMeth[DMR_fem_percMeth$ID %in% ac.id,]
DMR_fem_percMeth.c <- DMR_fem_percMeth[DMR_fem_percMeth$ID %in% c.id,]
DMR_mal_percMeth.ac <- DMR_mal_percMeth[DMR_mal_percMeth$ID %in% ac.id,]
DMR_mal_percMeth.c <- DMR_mal_percMeth[DMR_mal_percMeth$ID %in% c.id,]

# check for near 0 variance and remove
if (length(nearZeroVar(DMS_fem_percMeth.ac, freqCut = 15, uniqueCut = 20)) > 0) {
  DMS_fem_percMeth.ac <- DMS_fem_percMeth.ac[, -nearZeroVar(DMS_fem_percMeth.ac, freqCut = 15, uniqueCut = 20)] 
}

if (length(nearZeroVar(DMS_fem_percMeth.c, freqCut = 15, uniqueCut = 20)) > 0) {
  DMS_fem_percMeth.c <- DMS_fem_percMeth.c[, -nearZeroVar(DMS_fem_percMeth.c, freqCut = 15, uniqueCut = 20)] 
}

if (length(nearZeroVar(DMS_mal_percMeth.ac, freqCut = 15, uniqueCut = 20)) > 0) {
  DMS_mal_percMeth.ac <- DMS_mal_percMeth.ac[, -nearZeroVar(DMS_mal_percMeth.ac, freqCut = 15, uniqueCut = 20)] 
}

if (length(nearZeroVar(DMS_mal_percMeth.c, freqCut = 15, uniqueCut = 20)) > 0) {
  DMS_mal_percMeth.c <- DMS_mal_percMeth.c[, -nearZeroVar(DMS_mal_percMeth.c, freqCut = 15, uniqueCut = 20)] 
}


if (length(nearZeroVar(DMR_fem_percMeth.ac, freqCut = 15, uniqueCut = 20)) > 0) {
  DMR_fem_percMeth.ac <- DMR_fem_percMeth.ac[, -nearZeroVar(DMR_fem_percMeth.ac, freqCut = 15, uniqueCut = 20)] 
}

if (length(nearZeroVar(DMR_fem_percMeth.c, freqCut = 15, uniqueCut = 20)) > 0) {
  DMR_fem_percMeth.c <- DMR_fem_percMeth.c[, -nearZeroVar(DMR_fem_percMeth.c, freqCut = 15, uniqueCut = 20)] 
}

if (length(nearZeroVar(DMR_mal_percMeth.ac, freqCut = 15, uniqueCut = 20)) > 0) {
  DMR_mal_percMeth.ac <- DMR_mal_percMeth.ac[, -nearZeroVar(DMR_mal_percMeth.ac, freqCut = 15, uniqueCut = 20)] 
}

if (length(nearZeroVar(DMR_mal_percMeth.c, freqCut = 15, uniqueCut = 20)) > 0) {
  DMR_mal_percMeth.c <- DMR_mal_percMeth.c[, -nearZeroVar(DMR_mal_percMeth.c, freqCut = 15, uniqueCut = 20)] 
}

#run
fem.shoal.DMS.ac <- enet.caret.function.v2(shoal.data.ac, 
                                           DMS_fem_percMeth.ac, 
                                           "total_diff_s", "gaussian")
mal.shoal.DMS.ac <- enet.caret.function.v2(shoal.data.ac, 
                                           DMS_mal_percMeth.ac, 
                                           "total_diff_s","gaussian")

fem.shoal.DMR.ac <- enet.caret.function.v2(shoal.data.ac, 
                                           DMR_fem_percMeth.ac, 
                                           "total_diff_s", "gaussian")
mal.shoal.DMR.ac <- enet.caret.function.v2(shoal.data[shoal.data$ID %in% ac.id,], 
                                           DMR_mal_percMeth.ac, 
                                           "total_diff_s","gaussian")

fem.shoal.DMS.c <- enet.caret.function.v2(shoal.data.c, 
                                          DMS_fem_percMeth.c, 
                                          "total_diff_s", "gaussian")
mal.shoal.DMS.c <- enet.caret.function.v2(shoal.data.c, 
                                          DMS_mal_percMeth.c, 
                                          "total_diff_s","gaussian")

fem.shoal.DMR.c <- enet.caret.function.v2(shoal.data.c, 
                                          DMR_fem_percMeth.c, 
                                          "total_diff_s", "gaussian")
mal.shoal.DMR.c <- enet.caret.function.v2(shoal.data.c, 
                                          DMR_mal_percMeth.c, 
                                          "total_diff_s","gaussian")

## RESUME here so read back in data ##
fem.shoal.DMS.old <- readRDS("./dev_exp/data/enet_res/fem_shoal_DMS.RDS")
mal.shoal.DMS.old <- readRDS("./dev_exp/data/enet_res/mal_shoal_DMS.RDS")
fem.shoal.DMR.old <- readRDS("./dev_exp/data/enet_res/fem_shoal_DMR.RDS")
mal.shoal.DMR.old <- readRDS("./dev_exp/data/enet_res/mal_shoal_DMR.RDS")

#get top coefficients 
fem.shoal.DMS.top.coef <- get.top.coefficients(fem.shoal.DMS.old$model, top_n = 10)
mal.shoal.DMS.top.coef <- get.top.coefficients(mal.shoal.DMS.old$model, top_n = 10)

fem.shoal.DMR.top.coef <- get.top.coefficients(fem.shoal.DMR.old$model, top_n = 10)
mal.shoal.DMR.top.coef <- get.top.coefficients(mal.shoal.DMR.old$model, top_n = 10)

fem.shoal.DMS.top.coef.ac <- get.top.coefficients(fem.shoal.DMS.ac, top_n = 10)
mal.shoal.DMS.top.coef.ac <- get.top.coefficients(mal.shoal.DMS.ac, top_n = 10)

fem.shoal.DMR.top.coef.ac <- get.top.coefficients(fem.shoal.DMR.ac, top_n = 10)
mal.shoal.DMR.top.coef.ac <- get.top.coefficients(mal.shoal.DMR.ac, top_n = 10)

fem.shoal.DMS.top.coef.c <- get.top.coefficients(fem.shoal.DMS.c, top_n = 10)
mal.shoal.DMS.top.coef.c <- get.top.coefficients(mal.shoal.DMS.c, top_n = 10)

fem.shoal.DMR.top.coef.c <- get.top.coefficients(fem.shoal.DMR.c, top_n = 10)
mal.shoal.DMR.top.coef.c <- get.top.coefficients(mal.shoal.DMR.c, top_n = 10)

#run linear models
#DMSs
fem.shoal.DMS.lm <- run.lm.func(DMS_fem_percMeth, shoal.data, fem.shoal.DMS.top.coef$variable,"total_diff_s")
mal.shoal.DMS.lm <- run.lm.func(DMS_mal_percMeth, shoal.data, mal.shoal.DMS.top.coef$variable,"total_diff_s")

#DMRs
fem.shoal.DMR.lm <- run.lm.func(DMR_fem_percMeth, shoal.data, fem.shoal.DMR.top.coef$variable,"total_diff_s")
mal.shoal.DMR.lm <- run.lm.func(DMR_mal_percMeth, shoal.data, mal.shoal.DMR.top.coef$variable,"total_diff_s")

#DHARMa validation for all models 
fem.shoal.DMS.val.plot <- val.func(fem.shoal.DMS.lm)
mal.shoal.DMS.val.plot <- val.func(mal.shoal.DMS.lm)

fem.shoal.DMR.val.plot <- val.func(fem.shoal.DMR.lm)
mal.shoal.DMR.val.plot <- val.func(mal.shoal.DMR.lm)

#plot significant DMRs 
#merge datasets
fem.shoal.DMS.data <- merge(DMS_fem_percMeth, shoal.data, by = "ID")
fem.shoal.DMS.data$NC_024347.1_3570681 <- scale(fem.shoal.DMS.data$NC_024347.1_3570681)
fem.shoal.DMS.data$NC_024352.1_6869839 <- scale(fem.shoal.DMS.data$NC_024352.1_6869839)

mal.shoal.DMS.data <- merge(DMS_mal_percMeth, shoal.data, by = "ID")
mal.shoal.DMS.data$NC_024344.1_21712760 <- scale(mal.shoal.DMS.data$NC_024344.1_21712760)

mal.shoal.DMR.data <- merge(DMR_mal_percMeth, shoal.data, by = "ID")
mal.shoal.DMR.data$NC_024334.1_18706701 <- scale(mal.shoal.DMR.data$NC_024334.1_18706701)
mal.shoal.DMR.data$NC_024344.1_7298901 <- scale(mal.shoal.DMR.data$NC_024344.1_7298901)
mal.shoal.DMR.data$NC_024347.1_12981601 <- scale(mal.shoal.DMR.data$NC_024347.1_12981601)

fem.shoal.DMR.data <- merge(DMR_fem_percMeth, shoal.data, by = "ID")
fem.shoal.DMR.data$NC_024341.1_18706501 <- scale(fem.shoal.DMR.data$NC_024341.1_18706501)

#predict variables  
pred.mal.shoal.NC_024334.1_18706701 <- ggpredict(mal.shoal.DMR.lm, terms = c("NC_024334.1_18706701"))
pred.mal.shoal.NC_024344.1_7298901 <- ggpredict(mal.shoal.DMR.lm, terms = c("NC_024344.1_7298901"))
pred.mal.shoal.NC_024347.1_12981601 <- ggpredict(mal.shoal.DMR.lm, terms = c("NC_024347.1_12981601"))

pred.fem.shoal.NC_024341.1_18706501 <- ggpredict(fem.shoal.DMR.lm, terms = c("NC_024341.1_18706501"))

pred.fem.shoal.NC_024347.1_3570681 <- ggpredict(fem.shoal.DMS.lm, terms = c("NC_024347.1_3570681"))
pred.fem.shoal.NC_024352.1_6869839 <- ggpredict(fem.shoal.DMS.lm, terms = c("NC_024352.1_6869839"))

pred.mal.shoal.NC_024344.1_21712760 <- ggpredict(mal.shoal.DMS.lm, terms = c("NC_024344.1_21712760"))

#plot
NC_024341.1_18706501.shoal.plot <- ggplot(pred.fem.shoal.NC_024341.1_18706501) + 
  geom_point(data = fem.shoal.DMR.data, aes(x = NC_024341.1_18706501, y = total_diff_s)) + 
  xlab("NC_024341.1_18706501 Meth.") + ylab("Preference for shoaling (s)") +
  theme_bw() + geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error),
                           fill = "skyblue", alpha = 0.5)  +
  geom_line(aes(x = x, y = predicted), col = "black") + 
  theme(axis.text = element_text(size=13), axis.title = element_text(size=15)) 
NC_024341.1_18706501.shoal.plot

NC_024334.1_18706701.shoal.plot <- ggplot(pred.mal.shoal.NC_024334.1_18706701) + 
  geom_point(data = mal.shoal.DMR.data, aes(x = NC_024334.1_18706701, y = total_diff_s)) + 
  xlab("NC_024334.1_18706701 Meth.") + ylab("Preference for shoaling (s)") +
  theme_bw() + geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error),
                           fill = "skyblue", alpha = 0.5)  +
  geom_line(aes(x = x, y = predicted), col = "black") + 
  theme(axis.text = element_text(size=13), axis.title = element_text(size=15)) 

NC_024334.1_18706701.shoal.plot

NC_024344.1_72989011.shoal.plot <- ggplot(pred.mal.shoal.NC_024344.1_7298901) + 
  geom_point(data = mal.shoal.DMR.data, aes(x = NC_024344.1_7298901, y = total_diff_s)) + 
  xlab("NC_024344.1_7298901 Meth.") + ylab("Preference for shoaling (s)") +
  theme_bw() + geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error),
                           fill = "skyblue", alpha = 0.5)  +
  geom_line(aes(x = x, y = predicted), col = "black") + 
  theme(axis.text = element_text(size=13), axis.title = element_text(size=15)) 
NC_024344.1_72989011.shoal.plot

NC_024347.1_12981601.shoal.plot <- ggplot(pred.mal.shoal.NC_024347.1_12981601) + 
  geom_point(data = mal.shoal.DMR.data, aes(x = NC_024347.1_12981601, y = total_diff_s)) + 
  xlab("NC_024347.1_12981601 Meth.") + ylab("Preference for shoaling (s)") +
  theme_bw() + geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error),
                           fill = "skyblue", alpha = 0.5)  +
  geom_line(aes(x = x, y = predicted), col = "black") + 
  theme(axis.text = element_text(size=13), axis.title = element_text(size=15)) 
NC_024347.1_12981601.shoal.plot



NC_024347.1_3570681.shoal.plot <- ggplot(pred.fem.shoal.NC_024347.1_3570681) + 
  geom_point(data = fem.shoal.DMS.data, aes(x = NC_024347.1_3570681, y = total_diff_s)) + 
  xlab("NC_024347.1_3570681 Meth.") + ylab("Preference for shoal (s)") +
  theme_bw() + geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error),
                           fill = "skyblue", alpha = 0.5)  +
  geom_line(aes(x = x, y = predicted), col = "black") + 
  theme(axis.text = element_text(size=13), axis.title = element_text(size=15)) 
NC_024347.1_3570681.shoal.plot

NC_024352.1_6869839.shoal.plot <- ggplot(pred.fem.shoal.NC_024352.1_6869839) + 
  geom_point(data = fem.shoal.DMS.data, aes(x = NC_024352.1_6869839, y = total_diff_s)) + 
  xlab("NC_024352.1_6869839 Meth.") + ylab("Preference for shoal (s)") +
  theme_bw() + geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error),
                           fill = "skyblue", alpha = 0.5)  +
  geom_line(aes(x = x, y = predicted), col = "black") + 
  theme(axis.text = element_text(size=13), axis.title = element_text(size=15)) 
NC_024352.1_6869839.shoal.plot

NC_024344.1_21712760.shoal.plot <- ggplot(pred.mal.shoal.NC_024344.1_21712760) + 
  geom_point(data = mal.shoal.DMS.data, aes(x = NC_024344.1_21712760, y = total_diff_s)) + 
  xlab("NC_024344.1_21712760 Meth.") + ylab("Preference for shoal (s)") +
  theme_bw() + geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error),
                           fill = "skyblue", alpha = 0.5)  +
  geom_line(aes(x = x, y = predicted), col = "black") + 
  theme(axis.text = element_text(size=13), axis.title = element_text(size=15)) 
NC_024344.1_21712760.shoal.plot

#make panels
panel.dmr <- ggarrange(NC_024341.1_18706501.shoal.plot, NC_024334.1_18706701.shoal.plot, NC_024344.1_72989011.shoal.plot, NC_024347.1_12981601.shoal.plot,
                   labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
panel.dmr

tiff("lm_DMR_shoal_panel.tiff", units="in", width = 8, height = 8, res = 600)
panel.dmr
dev.off()

panel.dms <- ggarrange(NC_024347.1_3570681.shoal.plot, NC_024352.1_6869839.shoal.plot, NC_024344.1_21712760.shoal.plot,
                       labels = c("A", "B", "C"), ncol = 2, nrow = 2)

panel.dms

tiff("lm_DMS_shoal_panel.tiff", units="in", width = 8, height = 8, res = 600)
panel.dms
dev.off()

#save results 

## Investigate repeated DMRs/DMSs ##
#check for repeats in top DMRs/DMSs
fem.DMS <- c(fem.act.DMS.top.coef$variable,fem.bold.DMS.top.coef$variable,fem.exp.DMS.top.coef$variable,fem.shoal.DMS.top.coef$variable)
fem.DMR <- c(fem.act.DMR.top.coef$variable,fem.bold.DMR.top.coef$variable,fem.exp.DMR.top.coef$variable,fem.shoal.DMR.top.coef$variable)

mal.DMS <- c(mal.act.DMS.top.coef$variable,mal.bold.DMS.top.coef$variable,mal.shoal.DMS.top.coef$variable)
mal.DMR <- c(mal.act.DMR.top.coef$variable,mal.bold.DMR.top.coef$variable,mal.shoal.DMR.top.coef$variable)

table(fem.DMS)
table(fem.DMR)

table(mal.DMS)
table(mal.DMR)

#check for overlap between sexes
intersect(fem.DMS, mal.DMS)
intersect(fem.DMR, mal.DMR)

#find locations of repeated DMSs DMRs #
#read in repeated datasets 
ref.anno <- readTranscriptFeatures("./dev_exp/data/annotation/guppy.sorted.bed", unique.prom = FALSE, up.flank = 1500, down.flank = 500)

repeated.fem.DMS <- data.frame(chr = c("NC_024334.1", "NC_024340.1", "NC_024340.1", "NC_024345.1", "NC_024346.1", "NC_024348.1"),
                               start = c(17007057, 3038245, 31459107, 3369212, 27624493, 6976471),
                               end = c(17007057, 3038245, 31459107, 3369212, 27624493, 6976471))

repeated.fem.DMR <- data.frame(chr = c("NC_024335.1", "NC_024335.1", "NC_024336.1", "NC_024340.1","NC_024341.1","NC_024344.1",
                                       "NC_024345.1", "NC_024345.1","NC_024348.1"),
                               start = c(24281001,2977701,19209201,14972801,18706501,9921201,10692801,4148501,1546501),
                               end = c(24281100,2977800,19209300,14972900,18706600,9921300,10692900,4148600,1546600))

repeated.mal.DMS <- data.frame(chr = c("NC_024346.1", "NC_024353.1"),
                               start = c(28697188,13933874),
                               end = c(28697188,13933874))

repeated.mal.DMR <- data.frame(chr = c("NC_024333.1","NC_024334.1","NC_024338.1","NC_024339.1","NC_024339.1","NC_024339.1",
                                       "NC_024348.1","NC_024352.1"),
                               start = c(33015601,18706701,10300001,10231201,14719501,9849301,5438801,7403301),
                               end = c(33015700,18706800,10300100,10231300,14719600,9849400,5438900,7403400))

#rename chromosomes in each GRanges
repeated.fem.DMS.gr <- renameChr(repeated.fem.DMS)
repeated.fem.DMR.gr <- renameChr(repeated.fem.DMR)

repeated.mal.DMS.gr <- renameChr(repeated.mal.DMS)
repeated.mal.DMR.gr <- renameChr(repeated.mal.DMR)

## Annotate ##
#annotate with gene parts 
repeated.DMS.fem.anno <- annotateWithGeneParts(repeated.fem.DMS.gr, ref.anno)
repeated.DMR.fem.anno <- annotateWithGeneParts(repeated.fem.DMR.gr, ref.anno)

repeated.DMS.mal.anno <- annotateWithGeneParts(repeated.mal.DMS.gr, ref.anno)
repeated.DMR.mal.anno <- annotateWithGeneParts(repeated.mal.DMR.gr, ref.anno)

#get nearest TSS for DMS and DMRs
repeated.DMS.fem.tss <- getAssociationWithTSS(repeated.DMS.fem.anno)
repeated.DMR.fem.tss <- getAssociationWithTSS(repeated.DMR.fem.anno)

repeated.DMS.mal.tss <- getAssociationWithTSS(repeated.DMS.mal.anno)
repeated.DMR.mal.tss <- getAssociationWithTSS(repeated.DMR.mal.anno)

#get proper biomart
bm <- useEnsembl(biomart = "ensembl")
bm <- useDataset("preticulata_gene_ensembl", mart = bm)

#grab go terms
go_list <- getBM(mart = bm, attributes = c('ensembl_transcript_id','external_gene_name','go_id', "name_1006",
                                           "namespace_1003", "go_linkage_type"))
#prep go terms for gene mappings
go_list_data <- go_list[,c(3,6,1)] #limit to needed data
go_list_data <- go_list_data[go_list_data$go_id != '',] #remove empty go terms

