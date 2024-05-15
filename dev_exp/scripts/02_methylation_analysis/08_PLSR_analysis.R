#######################################################
### Goal: Plot venn diagram and heatmaps and do cluster analysis
### Author: Janay Fox
### R script
#######################################################

## Set up ##
#install packages
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("methylKit")
#install.packages("pls")

library(ggplot2)
library(dplyr)
library(methylKit)
library(pls)
library(lme4)
library(car)
library(plsRglm)

#read behaviour data 
of.data <- read.csv("./dev_exp/data/clean/clean_openfield.csv")
sh.data <- read.csv("./dev_exp/data/clean/clean_shoaling.csv")

#change ID name to match 
of.data$ID <- paste0("D", of.data$ID)
sh.data$ID <- paste0("D", sh.data$ID)

#calculate time differences for shoal vs empty container
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
  
  return(t.df)
}

#Prepare data for PLSR ##
#get percent methylation 
DMS_all_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMS_res/DMSdiffmeth_all_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMS_res/DMSmeth_all_5X.RDS")
DMS_all_percMeth <- prep_data(DMS_all_percMeth)

#get perc meth for DMRs 
DMR_all_percMeth <- get_percMeth_matrix("./dev_exp/data/methylkit_res/DMR_res/DMRdiffmeth_all_5X.RDS",
                                        "./dev_exp/data/methylkit_res/DMR_res/DMR_tile_meth_all_5X.RDS")
DMR_all_percMeth <- prep_data(DMR_all_percMeth)

# #get percent methylation for all CpGs
# #read data files
# CpG_meth_all <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSmeth_all_5X.RDS")
# CpG_meth_fem <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSmeth_fem_5X.RDS")
# CpG_meth_mal <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSmeth_mal_5X.RDS")
# 
# #calculate percent methylation 
# CpG_all_percMeth <- as.data.frame(percMethylation(CpG_meth_all))
# CpG_fem_percMeth <- as.data.frame(percMethylation(CpG_meth_fem))
# CpG_mal_percMeth <- as.data.frame(percMethylation(CpG_meth_mal))
# 
# #add site names
# CpG_meth_all <- nameSite(CpG_meth_all)
# CpG_meth_fem <- nameSite(CpG_meth_fem)
# CpG_meth_mal <- nameSite(CpG_meth_mal)
# 
# #add rownames to perc meth matrix  
# rownames(CpG_all_percMeth) <- CpG_meth_all$site_name
# rownames(CpG_fem_percMeth) <- CpG_meth_fem$site_name
# rownames(CpG_mal_percMeth) <- CpG_meth_mal$site_name
# 
# #remove extra file 
# rm(CpG_meth_all)
# rm(CpG_meth_fem)
# rm(CpG_meth_mal)
# 
# CpG_all_percMeth <- prep_data(CpG_all_percMeth)
# CpG_fem_percMeth <- prep_data(CpG_fem_percMeth)
# CpG_mal_percMeth <- prep_data(CpG_mal_percMeth)

#extract variables for analysis and scale of data 
of.data <- of.data[,c(1,5,11)]

sh.data.diff <- sh.data[,c(1,34)]

sh.data.cue.dat <- sh.data[,c(1:4,34)]

#prep for binding 
row.names(of.data) <- of.data$ID
of.data <- of.data[, -1]
of.data <- scale(of.data)

row.names(sh.data.diff) <- sh.data.diff$ID
sh.data.diff <- sh.data.diff[, -1, drop=FALSE]

row.names(sh.data.cue.dat) <- sh.data.cue.dat$ID
sh.data.cue.dat <- sh.data.cue.dat[, -1, drop=FALSE]

#convert to dataframe for binding
DMS_all_percMeth <- as.data.frame(DMS_all_percMeth)
DMR_all_percMeth <- as.data.frame(DMR_all_percMeth)

#bind data 
DMS_of_all_data <- merge(DMS_all_percMeth, of.data, by = "row.names", all = FALSE)
rownames(DMS_of_all_data) <- DMS_of_all_data[,1]
DMS_of_all_data[,1] <- NULL

DMS_sh_all_data <- merge(DMS_all_percMeth, sh.data.diff, by = "row.names", all = FALSE)
rownames(DMS_sh_all_data) <- DMS_sh_all_data[,1]
DMS_sh_all_data[,1] <- NULL

DMR_of_all_data <- merge(DMR_all_percMeth, of.data, by = "row.names", all = FALSE)
rownames(DMR_of_all_data) <- DMR_of_all_data[,1]
DMR_of_all_data[,1] <- NULL

DMR_sh_all_data <- merge(DMR_all_percMeth, sh.data.diff, by = "row.names", all = FALSE)
rownames(DMR_sh_all_data) <- DMR_sh_all_data[,1]
DMR_sh_all_data[,1] <- NULL

DMR_all_data <- merge(DMR_of_all_data, sh.data.cue.dat, by = "row.names", all = FALSE)
rownames(DMR_all_data) <- DMR_all_data[,1]
DMR_all_data[,1] <- NULL

## Try correlation analysis ## 
DMR.sh.cor <- as.data.frame(cor(DMR_sh_all_data))
DMR.of.cor <- as.data.frame(cor(DMR_of_all_data))

#subset based on corr  
DMR.sh.cor.sub <- DMR.sh.cor[DMR.sh.cor$total_diff_s > 0.1 | DMR.sh.cor$total_diff_s < -0.1 ,]
DMR.of.dist.cor.sub <- DMR.of.cor[DMR.of.cor$dist_cm > 0.1 | DMR.of.cor$dist_cm < -0.1 ,]
DMR.of.shelt.cor.sub <- DMR.of.cor[DMR.of.cor$time_shelter_s > 0.1 | DMR.of.cor$time_shelter_s < -0.1 ,]

#make lists of DMRs to run lmms on for each behaviour
sh_DMR_list <- rownames(DMR.sh.cor.sub)
dist_DMR_list <- rownames(DMR.of.dist.cor.sub)
shelt_DMR_list <- rownames(DMR.of.shelt.cor.sub)

#remove behavioural measurements from list
sh_DMR_list <- head(sh_DMR_list, -1)
dist_DMR_list <- head(dist_DMR_list, -2)
shelt_DMR_list <- head(shelt_DMR_list, -2)

#check correlations between DMRs 
cor(DMR_sh_all_data[,sh_DMR_list, drop = FALSE])

## LMM analysis ##
#function that runs lmms
model.function <- function(behav, DMR, data){
  
  data$DMR_asin <- asin(sqrt(data[[DMR]]/100))
  
  #run lmm
  lmm <- lmer(get(behav) ~ DMR_asin + sex + cue + (1|tank), REML = TRUE, data = data)

  print(summary(lmm))

  return(lmm)
}

#run lmms for shoaling
sh_lmm_list <- list()
sh_anova_res <- list()

for (DMR in sh_DMR_list) {
  #run model
  lmm <- model.function("total_diff_s", DMR, DMR_all_data)
  #save model
  sh_lmm_list[[DMR]] <- lmm
  #run anova
  anova_res <- Anova(lmm, type = 3)
  #save anova results
  sh_anova_res[[DMR]] <- anova_res
}

#correct p values
sh_fdr_corrected_list <- list()
for (DMR in names(sh_anova_res)) { 
  #extract the Anova result from the list
  anova_result <- sh_anova_res[[DMR]]
  
  #extract the uncorrected p-values
  uncorrected_pvalues <- anova_result$"Pr(>Chisq)"
  
  # Perform FDR correction using the Benjamini-Hochberg procedure
  fdr_corrected_pvalues <- p.adjust(uncorrected_pvalues, method = "BH")
  
  # Store the FDR-corrected p-values in the list
  sh_fdr_corrected_list[[DMR]] <- fdr_corrected_pvalues
}

#run lmms for distance
dist_lmm_list <- list()
dist_anova_res <- list()

for (DMR in dist_DMR_list) {
  #run model
  lmm <- model.function("dist_cm", DMR, DMR_all_data)
  #save model
  dist_lmm_list[[DMR]] <- lmm
  #run anova
  anova_res <- Anova(lmm, type = 3)
  #save anova results
  dist_anova_res[[DMR]] <- anova_res
}

#correct p values
dist_fdr_corrected_list <- list()
for (DMR in names(dist_anova_res)) { 
  #extract the Anova result from the list
  anova_result <- dist_anova_res[[DMR]]
  
  #extract the uncorrected p-values
  uncorrected_pvalues <- anova_result$"Pr(>Chisq)"
  
  # Perform FDR correction using the Benjamini-Hochberg procedure
  fdr_corrected_pvalues <- p.adjust(uncorrected_pvalues, method = "BH")
  
  # Store the FDR-corrected p-values in the list
  dist_fdr_corrected_list[[DMR]] <- fdr_corrected_pvalues
}

#run lmms for shelter
shelt_lmm_list <- list()
shelt_anova_res <- list()

for (DMR in shelt_DMR_list) {
  #run model
  lmm <- model.function("time_shelter_s", DMR, DMR_all_data)
  #save model
  shelt_lmm_list[[DMR]] <- lmm
  #run anova
  anova_res <- Anova(lmm, type = 3)
  #save anova results
  shelt_anova_res[[DMR]] <- anova_res
}

#correct p values
shelt_fdr_corrected_list <- list()
for (DMR in names(shelt_anova_res)) { 
  #extract the Anova result from the list
  anova_result <- shelt_anova_res[[DMR]]
  
  #extract the uncorrected p-values
  uncorrected_pvalues <- anova_result$"Pr(>Chisq)"
  
  # Perform FDR correction using the Benjamini-Hochberg procedure
  fdr_corrected_pvalues <- p.adjust(uncorrected_pvalues, method = "BH")
  
  # Store the FDR-corrected p-values in the list
  shelt_fdr_corrected_list[[DMR]] <- fdr_corrected_pvalues
}

## Fit PLSR models ##
#run transformation on all DMR/DMS data 
exclude_columns_of <- c("dist_cm", "time_shelter_s")
exclude_columns_sh <- "total_diff_s" 

#not sure if I want to use this or not?  
transformed_DMS_of <- DMS_of_all_data %>% mutate_at(vars(-exclude_columns_of), funs((asin(sqrt(./100)))))
transformed_DMR_of <- DMR_of_all_data %>% mutate_at(vars(-exclude_columns_of), funs((asin(sqrt(./100)))))
 
transformed_DMS_sh <- DMS_sh_all_data %>% mutate_at(vars(-exclude_columns_sh), funs((asin(sqrt(./100)))))
transformed_DMR_sh <- DMR_sh_all_data %>% mutate_at(vars(-exclude_columns_sh), funs((asin(sqrt(./100)))))

#run test models
#DMSs
DMS_dist_model_test <- plsr(dist_cm ~ ., data = dplyr::select(DMS_of_all_data, -time_shelter_s), 
                              ncomp = 20, validation = "CV", method = "simpls", scale = TRUE)

DMS_shelt_model_test <- plsr(time_shelter_s ~ ., data = dplyr::select(DMS_of_all_data, -dist_cm), 
                                ncomp = 20, validation = "CV", method = "simpls", scale = TRUE)

DMS_sh_model_test <- plsr(total_diff_s ~ ., data = DMS_sh_all_data, 
                              ncomp = 20, validation = "CV", method = "simpls", scale = TRUE)

#check for number of components 
summary(DMS_dist_model_test) #shows RMSEP and a bias corrected estimate
validationplot(DMS_dist_model_test) #plot 

summary(DMS_shelt_model_test)
validationplot(DMS_shelt_model_test) #plot 

summary(DMS_sh_model_test)
validationplot(DMS_sh_model_test) #plot 

#check plot 
plot(DMS_dist_all_model_test, asp = 1,  line = TRUE, ncomp = 5)
plot(DMS_dist_all_model_test, plottype = "scores", comps = 1:5)
plot(DMS_dist_all_model_test, plottype = "validation")

corrplot(DMS_dist_all_model_test, identify = TRUE)
loadingplot(DMS_dist_all_model_test)

#view results
co <- as.data.frame(coef(DMS_dist_all_model_test, ncomp = 2))
head(scores(DMS_dist_all_model_test, ncomp = 2))
head(loadings(DMS_dist_all_model_test, ncomp = 2))

R2(DMS_dist_all_model_test, estimate = "train")
R2(DMS_dist_all_model_test, estimate = "CV")

# #check results 
R2(DMS_of_all_model_test, ncomp = 5)
R2(DMS_sh_all_model_test, ncomp = 4)

r.squared
# #check models 

plot(DMS_sh_all_model_test, asp = 1,  line = TRUE, ncomp =4)

plot(DMR_sh_all_model_test, asp = 1,  line = TRUE, ncomp =2)

# #get regression coefficients (can do ncomp up to the num comps that I did in original model)
# DMS_of_fem_coef <- as.data.frame(coef(DMS_of_fem_model))
# DMS_sh_fem_coef <- coef(DMS_sh_fem_model)
# 
# DMS_of_fem_loadings <- as.data.frame(loadings(DMS_of_fem_model))
# DMS_sh_fem_loadings <- loadings(DMS_sh_fem_model)
# 
# DMS_of_fem_model$scores
# DMS_of_fem_model_test$validation$PRESS
# 
# 
# # # Extract Loadings
# # loadings_matrix <- loadings(DMS_of_fem_model)
# # 
# # # Sum of Squares of Loadings for Each Variable
# # ssq_loadings <- colSums(loadings_matrix^2)
# # 
# # # Variable Importance
# # vip_scores <- ssq_loadings / sum(ssq_loadings)
# 
# # 
# #predict test data
# DMS_of_fem_pls_pred <- as.data.frame(predict(DMS_of_fem_model, DMS_of_fem_test, ncomp = 4))
# 
# plot(DMS_of_fem_test$time_shelter_s, DMS_of_fem_pls_pred$`time_shelter_s.4 comps`, ylim=c(-11,2), xlim=c(-11,2))
# abline(0, 1, col="red")
# 
# cor.test(DMS_of_fem_test$time_shelter_s, DMS_of_fem_pls_pred$`time_shelter_s.4 comps`)
# 
# plot(DMS_sh_fem_model_test, ncomp = 4, asp = 1, line = TRUE)
# plot(DMS_sh_fem_model, plottype = "scores", comps = 1:4)
# 
# #this is variance of X so like how much variance in methylation
# explvar(DMS_sh_fem_model_test)
# 
# #predict
# predicted_data <- as.data.frame(predict(DMS_sh_fem_model))
# 
# RMSEP(DMS_sh_fem_model_test, newdata = predicted_data)
# 
# pred <- as.data.frame(predict(DMS_of_fem_model_1, ncomp = 10))
# 
# resp_dat <- DMS_of_fem_data[,c("dist_cm", "frozen_s", "time_shelter_s")]
# 
# test <- cor_coef <- cor(pred, resp_dat, method = "pearson")
# 
# p_values_matrix <- matrix(NA, nrow = ncol(resp_dat), ncol = ncol(pred))
# 
# varImp(DMS_sh_fem_model)
# 
# head(DMS_sh_fem_model$Yscores)
# head(DMS_sh_fem_model$model)
# #RMSE
# sqrt(mean((DMS_sh_fem_model$residuals)^2))
# 

DMR_sh_all_data$tank <- as.factor(DMR_sh_all_data$tank)
DMR_sh_all_data$cue <- as.factor(DMR_sh_all_data$cue)
DMR_sh_all_data$sex <- as.factor(DMR_sh_all_data$sex)

lm.test <- lmer(tight_diff_s ~ NC_024346.1_21149201 + cue + sex + (1|tank), data = DMR_sh_all_data)

plot(lm.test)

plot(tight_diff_s ~ NC_024346.1_21149201, data = DMR_sh_all_data)

summary(lm.test)

NC_024346.1_21149201
