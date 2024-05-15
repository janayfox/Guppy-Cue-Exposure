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
#install.packages("pls")

library(ggplot2)
library(dplyr)
library(methylKit)
library(lme4)
library(car)
library(plsRglm)
library(mixOmics)

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

#get data for all CpGs
#CpG_meth_all <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSmeth_all_5X.RDS")
CpG_meth_fem <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSmeth_fem_5X.RDS")
CpG_meth_mal <- readRDS("./dev_exp/data/methylkit_res/DMS_res/DMSmeth_mal_5X.RDS")

#calculate percent methylation
#perc_meth_CpG_all <- as.data.frame(percMethylation(CpG_meth_all))
CpG_fem_percMeth <- as.data.frame(percMethylation(CpG_meth_fem))
CpG_mal_percMeth <- as.data.frame(percMethylation(CpG_meth_mal))

#add site names
#CpG_meth_all <- nameSite(CpG_meth_all)
CpG_fem_percMeth <- nameSite(CpG_meth_fem)
CpG_mal_percMeth <- nameSite(CpG_meth_mal)

#add rownames to perc meth matrix
#rownames(perc_meth_CpG_all) <- CpG_meth_all$site_name
rownames(CpG_fem_percMeth) <- CpG_meth_fem$site_name
rownames(CpG_mal_percMeth) <- CpG_meth_mal$site_name

#remove extra file
#rm(CpG_meth_all)
rm(CpG_meth_fem)
rm(CpG_meth_mal)

#extract variables for analysis
of.data <- of.data[,c(1,5,11)]
sh.data.diff <- sh.data[,c(1,34)]
cue.dat <- sh.data[,c(1:4)]

#bind behavioural data 
behav.dat <- merge(of.data, sh.data.diff, by = "ID")
behav.dat <- merge(behav.dat, cue.dat, by = "ID")

#prep for binding with methylation data 
row.names(behav.dat) <- behav.dat$ID
behav.dat <- behav.dat[, -1]

#convert to dataframe for binding

# CpG_fem_percMeth <- as.data.frame(CpG_fem_percMeth)
DMS_fem_percMeth <- as.data.frame(DMS_fem_percMeth)
DMR_fem_percMeth <- as.data.frame(DMR_fem_percMeth)

# CpG_mal_percMeth <- as.data.frame(CpG_mal_percMeth)
DMS_mal_percMeth <- as.data.frame(DMS_mal_percMeth)
DMR_mal_percMeth <- as.data.frame(DMR_mal_percMeth)

#bind data 
DMS_fem_data <- merge(DMS_fem_percMeth, behav.dat, by = "row.names", all = FALSE)
rownames(DMS_fem_data) <- DMS_fem_data[,1]
DMS_fem_data[,1] <- NULL

DMR_fem_data <- merge(DMR_fem_percMeth, behav.dat, by = "row.names", fem = FALSE)
rownames(DMR_fem_data) <- DMR_fem_data[,1]
DMR_fem_data[,1] <- NULL

# CpG_mal_data <- merge(CpG_mal_percMeth, behav.dat, by = "row.names", all = FALSE)
# rownames(CpG_mal_data) <- CpG_mal_data[,1]
# CpG_mal_data[,1] <- NULL

DMS_mal_data <- merge(DMS_mal_percMeth, behav.dat, by = "row.names", mal = FALSE)
rownames(DMS_mal_data) <- DMS_mal_data[,1]
DMS_mal_data[,1] <- NULL

DMR_mal_data <- merge(DMR_mal_percMeth, behav.dat, by = "row.names", mal = FALSE)
rownames(DMR_mal_data) <- DMR_mal_data[,1]
DMR_mal_data[,1] <- NULL

## PLSR ## 
#run arcsin transformation DMR/DMS data 
exclude_columns_trans <- c("dist_cm", "time_shelter_s", "total_diff_s", "sex", "cue", "tank")

#asin_CpG_fem_data <- sapply(CpG_fem_percMeth, function(x) asin(sqrt(x/100)))
asin_DMS_fem_data <- DMS_fem_data %>% mutate_at(vars(-exclude_columns_trans), list(~ asin(sqrt(./100))))
asin_DMR_fem_data <- DMR_fem_data %>% mutate_at(vars(-exclude_columns_trans), list(~ asin(sqrt(./100))))

#asin_CpG_mal_data <- sapply(CpG_mal_percMeth, function(x) asin(sqrt(x/100)))
asin_DMS_mal_data <- DMS_mal_data %>% mutate_at(vars(-exclude_columns_trans), list(~ asin(sqrt(./100))))
asin_DMR_mal_data <- DMR_mal_data %>% mutate_at(vars(-exclude_columns_trans), list(~ asin(sqrt(./100))))

## mixOmics version ##
#DMRs
x_fem_DMR <- asin_DMR_fem_data[1:23]
x_mal_DMR <- asin_DMR_mal_data[1:206]

y_dist_fem <- DMS_fem_data[, "dist_cm"]
y_shelt_fem <- DMS_fem_data[, "time_shelter_s"]
y_shoal_fem <- DMS_fem_data[, "total_diff_s"]

y_dist_mal <- DMS_mal_data[, "dist_cm"]
y_shelt_mal <- DMS_mal_data[, "time_shelter_s"]
y_shoal_mal <- DMS_mal_data[, "total_diff_s"]

#check number of dimensions to use
#fem
fem.tune.dist.DMR <- pls(X = x_fem_DMR, Y = y_dist_fem, ncomp = 20, mode = 'regression') #make test plsr
fem.Q2.pls1.dist.DMR <- perf(fem.tune.dist.DMR, validation = 'Mfold', #calculate Q2 across components
                             folds = 10, nrepeat = 100, progressBar = TRUE)
plot(fem.Q2.pls1.dist.DMR, criterion = 'Q2') #visualize

fem.tune.shelt.DMR <- pls(X = x_fem_DMR, Y = y_shelt_fem, ncomp = 20, mode = 'regression') #make test plsr
fem.Q2.pls1.shelt.DMR <- perf(fem.tune.shelt.DMR, validation = 'Mfold', #calculate Q2 across components
                              folds = 10, nrepeat = 100, progressBar = TRUE)
plot(fem.Q2.pls1.shelt.DMR, criterion = 'Q2') #visualize

fem.tune.shoal.DMR <- pls(X = x_fem_DMR, Y = y_shoal_fem, ncomp = 20, mode = 'regression') #make test plsr
fem.Q2.pls1.shoal.DMR <- perf(fem.tune.shoal.DMR, validation = 'Mfold', #calculate Q2 across components
                              folds = 10, nrepeat = 100, progressBar = TRUE)
plot(fem.Q2.pls1.shoal.DMR, criterion = 'Q2') #visualize

#mal
mal.tune.dist.DMR <- pls(X = x_mal_DMR, Y = y_dist_mal, ncomp = 20, mode = 'regression') #make test plsr
mal.Q2.pls1.dist.DMR <- perf(mal.tune.dist.DMR, validation = 'Mfold', #calculate Q2 across components
                             folds = 10, nrepeat = 100, progressBar = TRUE)
plot(mal.Q2.pls1.dist.DMR, criterion = 'Q2') #visualize

mal.tune.shelt.DMR <- pls(X = x_mal_DMR, Y = y_shelt_mal, ncomp = 20, mode = 'regression') #make test plsr
mal.Q2.pls1.shelt.DMR <- perf(mal.tune.shelt.DMR, validation = 'Mfold', #calculate Q2 across components
                              folds = 10, nrepeat = 100, progressBar = TRUE)
plot(mal.Q2.pls1.shelt.DMR, criterion = 'Q2') #visualize

mal.tune.shoal.DMR <- pls(X = x_mal_DMR, Y = y_shoal_mal, ncomp = 20, mode = 'regression') #make test plsr
mal.Q2.pls1.shoal.DMR <- perf(mal.tune.shoal.DMR, validation = 'Mfold', #calculate Q2 across components
                              folds = 10, nrepeat = 100, progressBar = TRUE)
plot(mal.Q2.pls1.shoal.DMR, criterion = 'Q2') #visualize

#try splsr 
y_shelt_fem <- y_shelt_fem + 0.00001

list.keepX.fem <- c(2:20)
list.keepX.mal <- c(5:10, seq(15,50,5))

tune.spls.dist.fem <- tune.spls(X = x_fem_DMR, Y = y_dist_fem, ncomp = 2,
                                test.keepX = list.keepX.fem, validation = "Mfold",
                                folds = 10, nrepeat = 5, progressBar = TRUE,
                                measure = "MAE")
plot(tune.spls.dist.fem)
fem.dist.choice.ncomp <- tune.spls.dist.fem$choice.ncomp$ncomp
fem.dist.choice.keepX <- tune.spls.dist.fem$choice.keepX[1:fem.dist.choice.ncomp]
  
tune.spls.shelt.fem <- tune.spls(X = x_fem_DMR, Y = y_shelt_fem, ncomp = 2,
                                test.keepX = list.keepX.fem, validation = "Mfold",
                                folds = 10, nrepeat = 5, progressBar = TRUE,
                                measure = "MAE")
plot(tune.spls.shelt.fem)
fem.shelt.choice.ncomp <- tune.spls.shelt.fem$choice.ncomp$ncomp
fem.shelt.choice.keepX <- tune.spls.shelt.fem$choice.keepX[1:fem.dist.choice.ncomp]

tune.spls.shoal.fem <- tune.spls(X = x_fem_DMR, Y = y_shoal_fem, ncomp = 2,
                                test.keepX = list.keepX.fem, validation = "Mfold",
                                folds = 10, nrepeat = 5, progressBar = TRUE,
                                measure = "MAE")
plot(tune.spls.shoal.fem)
fem.shoal.choice.ncomp <- tune.spls.shoal.fem$choice.ncomp$ncomp
fem.shoal.choice.keepX <- tune.spls.shoal.fem$choice.keepX[1:fem.dist.choice.ncomp]

tune.spls.dist.mal <- tune.spls(X = x_mal_DMR, Y = y_dist_mal, ncomp = 2,
                                test.keepX = list.keepX.mal, validation = "Mfold",
                                folds = 10, nrepeat = 5, progressBar = TRUE,
                                measure = "MAE")
plot(tune.spls.dist.mal)
mal.dist.choice.ncomp <- tune.spls.dist.mal$choice.ncomp$ncomp
mal.dist.choice.keepX <- tune.spls.dist.mal$choice.keepX[1:fem.dist.choice.ncomp]

tune.spls.shelt.mal <- tune.spls(X = x_mal_DMR, Y = y_shelt_mal, ncomp = 2,
                                 test.keepX = list.keepX.mal, validation = "Mfold",
                                 folds = 10, nrepeat = 5, progressBar = TRUE,
                                 measure = "MAE")
plot(tune.spls.shelt.mal)
mal.shelt.choice.ncomp <- tune.spls.shelt.mal$choice.ncomp$ncomp
mal.shelt.choice.keepX <- tune.spls.shelt.mal$choice.keepX[1:fem.dist.choice.ncomp]

tune.spls.shoal.mal <- tune.spls(X = x_mal_DMR, Y = y_shoal_mal, ncomp = 2,
                                 test.keepX = list.keepX.mal, validation = "Mfold",
                                 folds = 10, nrepeat = 5, progressBar = TRUE,
                                 measure = "MAE")
plot(tune.spls.shoal.mal)
mal.shoal.choice.ncomp <- tune.spls.shoal.mal$choice.ncomp$ncomp
mal.shoal.choice.keepX <- tune.spls.shoal.mal$choice.keepX[1:fem.dist.choice.ncomp]

#run final spls moedls 
fem.spls.dist.DMR <- spls(X = x_fem_DMR, Y = y_dist_fem, ncomp = fem.dist.choice.ncomp, 
                          keepX = fem.dist.choice.keepX, mode = 'regression')
fem.spls.shelt.DMR <- spls(X = x_fem_DMR, Y = y_shelt_fem, ncomp = fem.shelt.choice.ncomp, 
                           keepX = fem.shelt.choice.keepX, mode = 'regression')
fem.spls.shoal.DMR <- spls(X = x_fem_DMR, Y = y_shoal_fem, ncomp = fem.shoal.choice.ncomp, 
                           keepX = fem.shoal.choice.keepX, mode = 'regression')

mal.spls.dist.DMR <- spls(X = x_mal_DMR, Y = y_dist_mal, ncomp = mal.dist.choice.ncomp, 
                          keepX = mal.dist.choice.keepX, mode = 'regression')
mal.spls.shelt.DMR <- spls(X = x_mal_DMR, Y = y_shelt_mal, ncomp = mal.shelt.choice.ncomp, 
                           keepX = mal.shelt.choice.keepX, mode = 'regression')
mal.spls.shoal.DMR <- spls(X = x_mal_DMR, Y = y_shoal_mal, ncomp = mal.shoal.choice.ncomp, 
                           keepX = mal.shoal.choice.keepX, mode = 'regression')

#run actual models
fem.pls.dist.DMR <- pls(X = x_fem_DMR, Y = y_dist_fem, ncomp = 1, mode = 'regression')
fem.pls.shelt.DMR <- pls(X = x_fem_DMR, Y = y_shelt_fem, ncomp = 1, mode = 'regression')
fem.pls.shoal.DMR <- pls(X = x_fem_DMR, Y = y_shoal_fem, ncomp = 1, mode = 'regression')

mal.pls.dist.DMR <- pls(X = x_mal_DMR, Y = y_dist_mal, ncomp = 1, mode = 'regression')
mal.pls.shelt.DMR <- pls(X = x_mal_DMR, Y = y_shelt_mal, ncomp = 1, mode = 'regression')
mal.pls.shoal.DMR <- pls(X = x_mal_DMR, Y = y_shoal_mal, ncomp = 1, mode = 'regression')

#calculate variation explained
fem.pls.dist.DMR$prop_expl_var$X
fem.pls.shelt.DMR$prop_expl_var$X
fem.pls.shoal.DMR$prop_expl_var$X

mal.pls.dist.DMR$prop_expl_var$X
mal.pls.shelt.DMR$prop_expl_var$X
mal.pls.shoal.DMR$prop_expl_var$X

#calculate correlation
cor.fem.dist <- cor(fem.pls.dist.DMR$variates$X, fem.pls.dist.DMR$variates$Y)[1,1]
cor.fem.shelt <- cor(fem.pls.shelt.DMR$variates$X, fem.pls.shelt.DMR$variates$Y)[1,1]
cor.fem.shoal <- cor(fem.pls.shoal.DMR$variates$X, fem.pls.shoal.DMR$variates$Y)[1,1]

cor.mal.dist <- cor(mal.pls.dist.DMR$variates$X, mal.pls.dist.DMR$variates$Y)[1,1]
cor.mal.shelt <- cor(mal.pls.shelt.DMR$variates$X, mal.pls.shelt.DMR$variates$Y)[1,1]
cor.mal.shoal <- cor(mal.pls.shoal.DMR$variates$X, mal.pls.shoal.DMR$variates$Y)[1,1]

cor.fem.dist.spls <- cor(fem.spls.dist.DMR$variates$X, fem.pls.dist.DMR$variates$Y)[1,1]
cor.fem.shelt.spls <- cor(fem.spls.shelt.DMR$variates$X, fem.pls.shelt.DMR$variates$Y)[1,1]
cor.fem.shoal.spls <- cor(fem.spls.shoal.DMR$variates$X, fem.pls.shoal.DMR$variates$Y)[1,1]

cor.mal.dist.spls <- cor(mal.spls.dist.DMR$variates$X, mal.pls.dist.DMR$variates$Y)[1,1]
cor.mal.shelt.spls <- cor(mal.spls.shelt.DMR$variates$X, mal.pls.shelt.DMR$variates$Y)[1,1]
cor.mal.shoal.spls <- cor(mal.spls.shoal.DMR$variates$X, mal.pls.shoal.DMR$variates$Y)[1,1]

#assss performace
fem.perf.dist <- perf(fem.pls.dist.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
fem.perf.shelt <- perf(fem.pls.shelt.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
fem.perf.shoal <- perf(fem.pls.shoal.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)

mal.perf.dist <- perf(mal.pls.dist.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
mal.perf.shelt <- perf(mal.pls.shelt.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
mal.perf.shoal <- perf(mal.pls.shoal.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)

fem.perf.dist.spls <- perf(fem.spls.dist.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
fem.perf.shelt.spls <- perf(fem.spls.shelt.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
fem.perf.shoal.spls <- perf(fem.spls.shoal.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)

mal.perf.dist.spls <- perf(mal.spls.dist.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
mal.perf.shelt.spls <- perf(mal.spls.shelt.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
mal.perf.shoal.spls <- perf(mal.spls.shoal.DMR, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)

#RMSEP
fem.perf.dist$measures$RMSEP
fem.perf.shelt$measures$RMSEP
fem.perf.shoal$measures$RMSEP

mal.perf.dist$measures$RMSEP
mal.perf.shelt$measures$RMSEP
mal.perf.shoal$measures$RMSEP

#R2 
fem.perf.dist$measures$R2
fem.perf.shelt$measures$R2
fem.perf.shoal$measures$R2

mal.perf.dist$measures$R2
mal.perf.shelt$measures$R2
mal.perf.shoal$measures$R2

fem.perf.dist.spls$measures$R2
fem.perf.shelt.spls$measures$R2
fem.perf.shoal.spls$measures$R2

mal.perf.dist.spls$measures$R2
mal.perf.shelt.spls$measures$R2
mal.perf.shoal.spls$measures$R2

#Assess significance 
## run permutations to assess significance ##
# Permutation test function
perm.func <- function(behavior_data, methylation_data, observed_correlation){
  n_permutations <- 100
  permuted_correlation <- numeric(n_permutations)
  
  for (i in 1:n_permutations) {
    # Permute behavior variables
    permuted_behavior <- sample(behavior_data)
    
    # Fit PLSR model with permuted behavior variables
    permuted_pls_model <- pls(Y = permuted_behavior, X = methylation_data, ncomp = 1, mode = 'regression')
    
    # Compute correlation with permuted behavior variables
    permuted_correlation[i] <- cor(permuted_pls_model$variates$X, permuted_pls_model$variates$Y)
  }
  
  # Calculate p-value
  p_value <- sum(permuted_correlation >= observed_correlation) / n_permutations
  print(p_value)
  return(p_value)
}

p_dist_fem <- perm.func(y_dist_fem, x_fem_DMR, cor.fem.dist)
p_shelt_fem <- perm.func(y_shelt_fem, x_fem_DMR, cor.fem.shelt)
p_shoal_fem <- perm.func(y_shoal_fem, x_fem_DMR, cor.fem.shoal)

p_dist_mal <- perm.func(y_dist_mal, x_mal_DMR, cor.mal.dist)
p_shelt_mal <- perm.func(y_shelt_mal, x_mal_DMR, cor.mal.shelt)
p_shoal_mal <- perm.func(y_shoal_mal, x_mal_DMR, cor.mal.shoal)

#spls perm function
spls.perm.func <- function(behavior_data, methylation_data, observed_correlation, ncomp_chc, nx){
  n_permutations <- 1000
  permuted_correlation <- numeric(n_permutations)
  
  for (i in 1:n_permutations) {
    # Permute behavior variables
    permuted_behavior <- sample(behavior_data)
    
    # Fit PLSR model with permuted behavior variables
    permuted_pls_model <- spls(Y = permuted_behavior, X = methylation_data, ncomp = ncomp_chc, 
                               keepX = nx, mode = 'regression')
    
    # Compute correlation with permuted behavior variables
    permuted_correlation[i] <- cor(permuted_pls_model$variates$X, permuted_pls_model$variates$Y)
  }
  
  # Calculate p-value
  p_value <- sum(permuted_correlation >= observed_correlation) / n_permutations
  print(p_value)
  return(p_value)
}

spls.perm.func(y_dist_fem, x_fem_DMR, cor.fem.dist.spls, fem.dist.choice.ncomp, fem.dist.choice.keepX)
spls.perm.func(y_shelt_fem, x_fem_DMR, cor.fem.shelt.spls, fem.shelt.choice.ncomp, fem.shelt.choice.keepX)
spls.perm.func(y_shoal_fem, x_fem_DMR, cor.fem.shoal.spls, fem.shoal.choice.ncomp, fem.shoal.choice.keepX)

spls.perm.func(y_dist_mal, x_mal_DMR, cor.mal.dist.spls, mal.dist.choice.ncomp, mal.dist.choice.keepX)
spls.perm.func(y_shelt_mal, x_mal_DMR, cor.mal.shelt.spls, mal.shelt.choice.ncomp, mal.shelt.choice.keepX)
spls.perm.func(y_shoal_mal, x_mal_DMR, cor.mal.shoal.spls, mal.shoal.choice.ncomp, mal.shoal.choice.keepX)

#look at VIP
fem.vip.pls.dist <- as.data.frame(vip(fem.spls.dist.DMR))
fem.vip.pls.shelt <- as.data.frame(vip(fem.spls.shelt.DMR))
fem.vip.pls.shoal <- as.data.frame(vip(fem.spls.shoal.DMR))

mal.vip.pls.dist <- as.data.frame(vip(mal.spls.dist.DMR))
mal.vip.pls.shelt <- as.data.frame(vip(mal.spls.shelt.DMR))
mal.vip.pls.shoal <- as.data.frame(vip(mal.spls.shoal.DMR))

#retain genes that have a VIP over 1
fem.vip.pls.dist <- rownames(subset(fem.vip.pls.dist, comp1 >= 1))
fem.vip.pls.shelt <- rownames(subset(fem.vip.pls.shelt, comp1 >= 1))
fem.vip.pls.shoal <- rownames(subset(fem.vip.pls.shoal, comp1 >= 1))

mal.vip.pls.dist <- rownames(subset(mal.vip.pls.dist, comp1 >= 1))
mal.vip.pls.shelt <- rownames(subset(mal.vip.pls.shelt, comp1 >= 1))
mal.vip.pls.shoal <- rownames(subset(mal.vip.pls.shoal, comp1 >= 1))

#DMSs
x_fem_DMS <- asin_DMS_fem_data[1:1269]
x_mal_DMS <- asin_DMS_mal_data[1:6930]

#check number of dimensions to use
#fem
fem.tune.dist.DMS <- pls(X = x_fem_DMS, Y = y_dist_fem, ncomp = 20, mode = 'regression') #make test plsr
fem.Q2.pls1.dist.DMS <- perf(fem.tune.dist.DMS, validation = 'Mfold', #calculate Q2 across components
                             folds = 10, nrepeat = 100, progressBar = TRUE)
plot(fem.Q2.pls1.dist.DMS, criterion = 'Q2') #visualize

fem.tune.shelt.DMS <- pls(X = x_fem_DMS, Y = y_shelt_fem, ncomp = 20, mode = 'regression') #make test plsr
fem.Q2.pls1.shelt.DMS <- perf(fem.tune.shelt.DMS, validation = 'Mfold', #calculate Q2 across components
                              folds = 10, nrepeat = 100, progressBar = TRUE)
plot(fem.Q2.pls1.shelt.DMS, criterion = 'Q2') #visualize

fem.tune.shoal.DMS <- pls(X = x_fem_DMS, Y = y_shoal_fem, ncomp = 20, mode = 'regression') #make test plsr
fem.Q2.pls1.shoal.DMS <- perf(fem.tune.shoal.DMS, validation = 'Mfold', #calculate Q2 across components
                              folds = 10, nrepeat = 100, progressBar = TRUE)
plot(fem.Q2.pls1.shoal.DMS, criterion = 'Q2') #visualize

#mal
mal.tune.dist.DMS <- pls(X = x_mal_DMS, Y = y_dist_mal, ncomp = 20, mode = 'regression') #make test plsr
mal.Q2.pls1.dist.DMS <- perf(mal.tune.dist.DMS, validation = 'Mfold', #calculate Q2 across components
                             folds = 10, nrepeat = 100, progressBar = TRUE)
plot(mal.Q2.pls1.dist.DMS, criterion = 'Q2') #visualize

mal.tune.shelt.DMS <- pls(X = x_mal_DMS, Y = y_shelt_mal, ncomp = 20, mode = 'regression') #make test plsr
mal.Q2.pls1.shelt.DMS <- perf(mal.tune.shelt.DMS, validation = 'Mfold', #calculate Q2 across components
                              folds = 10, nrepeat = 100, progressBar = TRUE)
plot(mal.Q2.pls1.shelt.DMS, criterion = 'Q2') #visualize

mal.tune.shoal.DMS <- pls(X = x_mal_DMS, Y = y_shoal_mal, ncomp = 20, mode = 'regression') #make test plsr
mal.Q2.pls1.shoal.DMS <- perf(mal.tune.shoal.DMS, validation = 'Mfold', #calculate Q2 across components
                              folds = 10, nrepeat = 100, progressBar = TRUE)
plot(mal.Q2.pls1.shoal.DMS, criterion = 'Q2') #visualize

#run actual models
fem.pls.dist.DMS <- pls(X = x_fem_DMS, Y = y_dist_fem, ncomp = 1, mode = 'regression')
fem.pls.shelt.DMS <- pls(X = x_fem_DMS, Y = y_shelt_fem, ncomp = 1, mode = 'regression')
fem.pls.shoal.DMS <- pls(X = x_fem_DMS, Y = y_shoal_fem, ncomp = 1, mode = 'regression')

mal.pls.dist.DMS <- pls(X = x_mal_DMS, Y = y_dist_mal, ncomp = 1, mode = 'regression')
mal.pls.shelt.DMS <- pls(X = x_mal_DMS, Y = y_shelt_mal, ncomp = 1, mode = 'regression')
mal.pls.shoal.DMS <- pls(X = x_mal_DMS, Y = y_shoal_mal, ncomp = 1, mode = 'regression')

#calculate variation explained
fem.pls.dist.DMS$prop_expl_var$X
fem.pls.shelt.DMS$prop_expl_var$X
fem.pls.shoal.DMS$prop_expl_var$X

mal.pls.dist.DMS$prop_expl_var$X
mal.pls.shelt.DMS$prop_expl_var$X
mal.pls.shoal.DMS$prop_expl_var$X

#assss performace
fem.perf.dist <- perf(fem.pls.dist.DMS, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
fem.perf.shelt <- perf(fem.pls.shelt.DMS, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
fem.perf.shoal <- perf(fem.pls.shoal.DMS, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)

mal.perf.dist <- perf(mal.pls.dist.DMS, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
mal.perf.shelt <- perf(mal.pls.shelt.DMS, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
mal.perf.shoal <- perf(mal.pls.shoal.DMS, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)

#RMSEP
fem.perf.dist$measures$RMSEP
fem.perf.shelt$measures$RMSEP
fem.perf.shoal$measures$RMSEP

mal.perf.dist$measures$RMSEP
mal.perf.shelt$measures$RMSEP
mal.perf.shoal$measures$RMSEP

#R2
fem.perf.dist$measures$R2
fem.perf.shelt$measures$R2
fem.perf.shoal$measures$R2

mal.perf.dist$measures$R2
mal.perf.shelt$measures$R2
mal.perf.shoal$measures$R2

#calculate correlation
cor.fem.dist.DMS <- cor(fem.pls.dist.DMS$variates$X, fem.pls.dist.DMS$variates$Y)[1,1]
cor.fem.shelt.DMS <- cor(fem.pls.shelt.DMS$variates$X, fem.pls.shelt.DMS$variates$Y)[1,1]
cor.fem.shoal.DMS <- cor(fem.pls.shoal.DMS$variates$X, fem.pls.shoal.DMS$variates$Y)[1,1]

cor.mal.dist.DMS <- cor(mal.pls.dist.DMS$variates$X, mal.pls.dist.DMS$variates$Y)[1,1]
cor.mal.shelt.DMS <- cor(mal.pls.shelt.DMS$variates$X, mal.pls.shelt.DMS$variates$Y)[1,1]
cor.mal.shoal.DMS <- cor(mal.pls.shoal.DMS$variates$X, mal.pls.shoal.DMS$variates$Y)[1,1]

#permutations 
p_dist_fem_DMS <- perm.func(y_dist_fem, x_fem_DMS, cor.fem.dist.DMS)
p_shelt_fem_DMS <- perm.func(y_shelt_fem, x_fem_DMS, cor.fem.shelt.DMS)
p_shoal_fem_DMS <- perm.func(y_shoal_fem, x_fem_DMS, cor.fem.shoal.DMS)

p_dist_mal_DMS <- perm.func(y_dist_mal, x_mal_DMS, cor.mal.dist.DMS)
p_shelt_mal_DMS <- perm.func(y_shelt_mal, x_mal_DMS, cor.mal.shelt.DMS)
p_shoal_mal_DMS <- perm.func(y_shoal_mal, x_mal_DMS, cor.mal.shoal.DMS)


#look at VIP
all.vip.pls.dist <- as.data.frame(vip(all.pls.dist.DMS))
all.vip.pls.shelt <- as.data.frame(vip(all.pls.shelt.DMS))
all.vip.pls.shoal <- as.data.frame(vip(all.pls.shoal.DMS))

fem.vip.pls.dist <- as.data.frame(vip(fem.pls.dist.DMS))
fem.vip.pls.shelt <- as.data.frame(vip(fem.pls.shelt.DMS))
fem.vip.pls.shoal <- as.data.frame(vip(fem.pls.shoal.DMS))

mal.vip.pls.dist <- as.data.frame(vip(mal.pls.dist.DMS))
mal.vip.pls.shelt <- as.data.frame(vip(mal.pls.shelt.DMS))
mal.vip.pls.shoal <- as.data.frame(vip(mal.pls.shoal.DMS))

#retain genes that have a VIP over 1
all.vip.pls.dist <- rownames(subset(all.vip.pls.dist, comp1 >= 1))
all.vip.pls.shelt <- rownames(subset(all.vip.pls.shelt, comp1 >= 1))
all.vip.pls.shoal <- rownames(subset(all.vip.pls.shoal, comp1 >= 1))

fem.vip.pls.dist <- rownames(subset(fem.vip.pls.dist, comp1 >= 1))
fem.vip.pls.shelt <- rownames(subset(fem.vip.pls.shelt, comp1 >= 1))
fem.vip.pls.shoal <- rownames(subset(fem.vip.pls.shoal, comp1 >= 1))

mal.vip.pls.dist <- rownames(subset(mal.vip.pls.dist, comp1 >= 1))
mal.vip.pls.shelt <- rownames(subset(mal.vip.pls.shelt, comp1 >= 1))
mal.vip.pls.shoal <- rownames(subset(mal.vip.pls.shoal, comp1 >= 1))

#CpGs
# x_fem_CpG <- asin_CpG_fem_data[1:23]
# x_mal_CpG <- asin_CpG_mal_data[1:206]
# 
# 
# #check number of dimensions to use
# #all
# all.tune.dist.CpG <- pls(X = x_all, Y = y_dist_all, ncomp = 20, mode = 'regression') #make test plsr
# all.Q2.pls1.dist.CpG <- perf(all.tune.dist.CpG, validation = 'Mfold', #calculate Q2 across components
#                              folds = 10, nrepeat = 100, progressBar = TRUE)
# plot(all.Q2.pls1.dist.CpG, criterion = 'Q2') #visualize
# all.choice.ncomp.dist <- 1 #select number of components
# 
# all.tune.shelt.CpG <- pls(X = x_all, Y = y_shelt_all, ncomp = 20, mode = 'regression') #make test plsr
# all.Q2.pls1.shelt.CpG <- perf(all.tune.shelt.CpG, validation = 'Mfold', #calculate Q2 across components
#                               folds = 10, nrepeat = 100, progressBar = TRUE)
# plot(all.Q2.pls1.shelt.CpG, criterion = 'Q2') #visualize
# all.choice.ncomp.shelt <- 1 #select number of components
# 
# all.tune.shoal.CpG <- pls(X = x_all, Y = y_shoal_all, ncomp = 20, mode = 'regression') #make test plsr
# all.Q2.pls1.shoal.CpG <- perf(all.tune.shoal.CpG, validation = 'Mfold', #calculate Q2 across components
#                               folds = 10, nrepeat = 100, progressBar = TRUE)
# plot(all.Q2.pls1.shoal.CpG, criterion = 'Q2') #visualize
# all.choice.ncomp.shoal <- 1 #select number of components
# 
# #fem 
# fem.tune.dist.CpG <- pls(X = x_fem, Y = y_dist_fem, ncomp = 20, mode = 'regression') #make test plsr
# fem.Q2.pls1.dist.CpG <- perf(fem.tune.dist.CpG, validation = 'Mfold', #calculate Q2 across components
#                              folds = 10, nrepeat = 100, progressBar = TRUE)
# plot(fem.Q2.pls1.dist.CpG, criterion = 'Q2') #visualize
# fem.choice.ncomp.dist <- 1 #select number of components
# 
# fem.tune.shelt.CpG <- pls(X = x_fem, Y = y_shelt_fem, ncomp = 20, mode = 'regression') #make test plsr
# fem.Q2.pls1.shelt.CpG <- perf(fem.tune.shelt.CpG, validation = 'Mfold', #calculate Q2 across components
#                               folds = 10, nrepeat = 100, progressBar = TRUE)
# plot(fem.Q2.pls1.shelt.CpG, criterion = 'Q2') #visualize
# fem.choice.ncomp.shelt <- 1 #select number of components
# 
# fem.tune.shoal.CpG <- pls(X = x_fem, Y = y_shoal_fem, ncomp = 20, mode = 'regression') #make test plsr
# fem.Q2.pls1.shoal.CpG <- perf(fem.tune.shoal.CpG, validation = 'Mfold', #calculate Q2 across components
#                               folds = 10, nrepeat = 100, progressBar = TRUE)
# plot(fem.Q2.pls1.shoal.CpG, criterion = 'Q2') #visualize
# fem.choice.ncomp.shoal <- 1 #select number of components
# 
# #mal
# mal.tune.dist.CpG <- pls(X = x_mal, Y = y_dist_mal, ncomp = 20, mode = 'regression') #make test plsr
# mal.Q2.pls1.dist.CpG <- perf(mal.tune.dist.CpG, validation = 'Mfold', #calculate Q2 across components
#                              folds = 10, nrepeat = 100, progressBar = TRUE)
# plot(mal.Q2.pls1.dist.CpG, criterion = 'Q2') #visualize
# mal.choice.ncomp.dist <- 1 #select number of components
# 
# mal.tune.shelt.CpG <- pls(X = x_mal, Y = y_shelt_mal, ncomp = 20, mode = 'regression') #make test plsr
# mal.Q2.pls1.shelt.CpG <- perf(mal.tune.shelt.CpG, validation = 'Mfold', #calculate Q2 across components
#                               folds = 10, nrepeat = 100, progressBar = TRUE)
# plot(mal.Q2.pls1.shelt.CpG, criterion = 'Q2') #visualize
# mal.choice.ncomp.shelt <- 1 #select number of components
# 
# mal.tune.shoal.CpG <- pls(X = x_mal, Y = y_shoal_mal, ncomp = 20, mode = 'regression') #make test plsr
# mal.Q2.pls1.shoal.CpG <- perf(mal.tune.shoal.CpG, validation = 'Mfold', #calculate Q2 across components
#                               folds = 10, nrepeat = 100, progressBar = TRUE)
# plot(mal.Q2.pls1.shoal.CpG, criterion = 'Q2') #visualize
# mal.choice.ncomp.shoal <- 1 #select number of components
# 
# #run actual models 
# all.pls.dist.CpG <- pls(X = x_all, Y = y_dist_all, ncomp = 1, mode = 'regression') 
# all.pls.shelt.CpG <- pls(X = x_all, Y = y_shelt_all, ncomp = 1, mode = 'regression') 
# all.pls.shoal.CpG <- pls(X = x_all, Y = y_shoal_all, ncomp = 1, mode = 'regression') 
# 
# fem.pls.dist.CpG <- pls(X = x_fem, Y = y_dist_fem, ncomp = 2, mode = 'regression') 
# fem.pls.shelt.CpG <- pls(X = x_fem, Y = y_shelt_fem, ncomp = 2, mode = 'regression') 
# fem.pls.shoal.CpG <- pls(X = x_fem, Y = y_shoal_fem, ncomp = 2, mode = 'regression') 
# 
# mal.pls.dist.CpG <- pls(X = x_mal, Y = y_dist_mal, ncomp = 2, mode = 'regression') 
# mal.pls.shelt.CpG <- pls(X = x_mal, Y = y_shelt_mal, ncomp = 2, mode = 'regression') 
# mal.pls.shoal.CpG <- pls(X = x_mal, Y = y_shoal_mal, ncomp = 2, mode = 'regression') 
# 
# #calculate variation explained 
# all.pls.dist.CpG$prop_expl_var$X
# all.pls.shelt.CpG$prop_expl_var$X
# all.pls.shoal.CpG$prop_expl_var$X
# 
# fem.pls.dist.CpG$prop_expl_var$X
# fem.pls.shelt.CpG$prop_expl_var$X
# fem.pls.shoal.CpG$prop_expl_var$X
# 
# mal.pls.dist.CpG$prop_expl_var$X
# mal.pls.shelt.CpG$prop_expl_var$X
# mal.pls.shoal.CpG$prop_expl_var$X
# 
# #calculate correlation
# cor(all.pls.dist.CpG$variates$X, all.pls.dist.CpG$variates$Y)
# cor(all.pls.shelt.CpG$variates$X, all.pls.shelt.CpG$variates$Y)
# cor(all.pls.shoal.CpG$variates$X, all.pls.shoal.CpG$variates$Y)
# 
# cor(fem.pls.dist.CpG$variates$X, fem.pls.dist.CpG$variates$Y)
# cor(fem.pls.shelt.CpG$variates$X, fem.pls.shelt.CpG$variates$Y)
# cor(fem.pls.shoal.CpG$variates$X, fem.pls.shoal.CpG$variates$Y)
# 
# cor(mal.pls.dist.CpG$variates$X, mal.pls.dist.CpG$variates$Y)
# cor(mal.pls.shelt.CpG$variates$X, mal.pls.shelt.CpG$variates$Y)
# cor(mal.pls.shoal.CpG$variates$X, mal.pls.shoal.CpG$variates$Y)
# 
# #assss performace 
# all.perf.dist <- perf(all.pls.dist.CpG, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
# all.perf.shelt <- perf(all.pls.shelt.CpG, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
# all.perf.shoal <- perf(all.pls.shoal.CpG, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
# 
# fem.perf.dist <- perf(fem.pls.dist.CpG, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
# fem.perf.shelt <- perf(fem.pls.shelt.CpG, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
# fem.perf.shoal <- perf(fem.pls.shoal.CpG, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
# 
# mal.perf.dist <- perf(mal.pls.dist.CpG, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
# mal.perf.shelt <- perf(mal.pls.shelt.CpG, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
# mal.perf.shoal <- perf(mal.pls.shoal.CpG, alidation = 'Mfold', folds = 10, nrepeat = 100, progressBar = TRUE)
# 
# all.perf.dist$measures$RMSEP
# all.perf.shelt$measures$RMSEP
# all.perf.shoal$measures$RMSEP
# 
# fem.perf.dist$measures$RMSEP
# fem.perf.shelt$measures$RMSEP
# fem.perf.shoal$measures$RMSEP
# 
# mal.perf.dist$measures$RMSEP
# mal.perf.shelt$measures$RMSEP
# mal.perf.shoal$measures$RMSEP
# 
# all.perf.dist$measures$R2
# all.perf.shelt$measures$R2
# all.perf.shoal$measures$R2
# 
# fem.perf.dist$measures$R2
# fem.perf.shelt$measures$R2
# fem.perf.shoal$measures$R2
# 
# mal.perf.dist$measures$R2
# mal.perf.shelt$measures$R2
# mal.perf.shoal$measures$R2
# 
# #worked better with 3 ??? 
# 
# 
# #look at VIP
# all.vip.pls.dist <- as.data.frame(vip(all.pls.dist.CpG))
# all.vip.pls.shelt <- as.data.frame(vip(all.pls.shelt.CpG))
# all.vip.pls.shoal <- as.data.frame(vip(all.pls.shoal.CpG))
# 
# fem.vip.pls.dist <- as.data.frame(vip(fem.pls.dist.CpG))
# fem.vip.pls.shelt <- as.data.frame(vip(fem.pls.shelt.CpG))
# fem.vip.pls.shoal <- as.data.frame(vip(fem.pls.shoal.CpG))
# 
# mal.vip.pls.dist <- as.data.frame(vip(mal.pls.dist.CpG))
# mal.vip.pls.shelt <- as.data.frame(vip(mal.pls.shelt.CpG))
# mal.vip.pls.shoal <- as.data.frame(vip(mal.pls.shoal.CpG))
# 
# #retain genes that have a VIP over 1
# all.vip.pls.dist <- rownames(subset(all.vip.pls.dist, comp1 >= 1))
# all.vip.pls.shelt <- rownames(subset(all.vip.pls.shelt, comp1 >= 1))
# all.vip.pls.shoal <- rownames(subset(all.vip.pls.shoal, comp1 >= 1))
# 
# fem.vip.pls.dist <- rownames(subset(fem.vip.pls.dist, comp1 >= 1))
# fem.vip.pls.shelt <- rownames(subset(fem.vip.pls.shelt, comp1 >= 1))
# fem.vip.pls.shoal <- rownames(subset(fem.vip.pls.shoal, comp1 >= 1))
# 
# mal.vip.pls.dist <- rownames(subset(mal.vip.pls.dist, comp1 >= 1))
# mal.vip.pls.shelt <- rownames(subset(mal.vip.pls.shelt, comp1 >= 1))
# mal.vip.pls.shoal <- rownames(subset(mal.vip.pls.shoal, comp1 >= 1))
# 
# ## run LMMs ##
# #run lmms for dist
# dist_lmm_list <- list()
# dist_anova_res <- list()
# 
# for (DMR in vip.pls.dist) {
#   #run model
#   #run lmm
#   lmm <- lmer(dist_cm ~ get(DMR) + sex + cue + (1|tank), REML = TRUE, data = asin_DMR_all_data)
#   
#   print(summary(lmm))
#   
#   #save model
#   dist_lmm_list[[DMR]] <- lmm
#   #run anova
#   anova_res <- Anova(lmm, type = 3)
#   #save anova results
#   dist_anova_res[[DMR]] <- anova_res
# }
# 
# #correct p values
# dist_fdr_corrected_list <- list()
# for (DMR in names(dist_anova_res)) { 
#   #extract the Anova result from the list
#   anova_result <- dist_anova_res[[DMR]]
#   
#   #extract the uncorrected p-values
#   uncorrected_pvalues <- anova_result$"Pr(>Chisq)"
#   
#   # Perform FDR correction using the Benjamini-Hochberg procedure
#   fdr_corrected_pvalues <- p.adjust(uncorrected_pvalues, method = "hommel")
#   
#   # Store the FDR-corrected p-values in the list
#   dist_fdr_corrected_list[[DMR]] <- fdr_corrected_pvalues
# }
# 
# #run lmms for dist
# dist_lmm_list <- list()
# dist_anova_res <- list()
# 
# for (DMR in vip.pls.dist) {
#   #run model
#   #run lmm
#   lmm <- lmer(dist_cm ~ get(DMR) + sex + cue + (1|tank), REML = TRUE, data = asin_DMR_all_data)
#   
#   print(summary(lmm))
#   
#   #save model
#   dist_lmm_list[[DMR]] <- lmm
#   #run anova
#   anova_res <- Anova(lmm, type = 3)
#   #save anova results
#   dist_anova_res[[DMR]] <- anova_res
# }
# 
# #correct p values
# dist_fdr_corrected_list <- list()
# for (DMR in names(dist_anova_res)) { 
#   #extract the Anova result from the list
#   anova_result <- dist_anova_res[[DMR]]
#   
#   #extract the uncorrected p-values
#   uncorrected_pvalues <- anova_result$"Pr(>Chisq)"
#   
#   # Perform FDR correction using the Benjamini-Hochberg procedure
#   fdr_corrected_pvalues <- p.adjust(uncorrected_pvalues, method = "hommel")
#   
#   # Store the FDR-corrected p-values in the list
#   dist_fdr_corrected_list[[DMR]] <- fdr_corrected_pvalues
# }
# 
# #run lmms for shelt
# shelt_lmm_list <- list()
# shelt_anova_res <- list()
# 
# for (DMR in vip.pls.shelt) {
#   #run model
#   #run lmm
#   lmm <- lmer(time_shelter_s ~ get(DMR) + sex + cue + (1|tank), REML = TRUE, data = asin_DMR_all_data)
#   
#   print(summary(lmm))
#   
#   #save model
#   shelt_lmm_list[[DMR]] <- lmm
#   #run anova
#   anova_res <- Anova(lmm, type = 3)
#   #save anova results
#   shelt_anova_res[[DMR]] <- anova_res
# }
# 
# #correct p values
# shelt_fdr_corrected_list <- list()
# for (DMR in names(shelt_anova_res)) { 
#   #extract the Anova result from the list
#   anova_result <- shelt_anova_res[[DMR]]
#   
#   #extract the uncorrected p-values
#   uncorrected_pvalues <- anova_result$"Pr(>Chisq)"
#   
#   # Perform FDR correction using the Benjamini-Hochberg procedure
#   fdr_corrected_pvalues <- p.adjust(uncorrected_pvalues, method = "hommel")
#   
#   # Store the FDR-corrected p-values in the list
#   shelt_fdr_corrected_list[[DMR]] <- fdr_corrected_pvalues
# }
# 
# #run lmms for shoal
# shoal_lmm_list <- list()
# shoal_anova_res <- list()
# 
# for (DMR in vip.pls.shoal) {
#   #run model
#   #run lmm
#   lmm <- lmer(total_diff_s ~ get(DMR) + sex + cue + (1|tank), REML = TRUE, data = asin_DMR_all_data)
#   
#   print(summary(lmm))
#   
#   #save model
#   shoal_lmm_list[[DMR]] <- lmm
#   #run anova
#   anova_res <- Anova(lmm, type = 3)
#   #save anova results
#   shoal_anova_res[[DMR]] <- anova_res
# }
# 
# #correct p values
# shoal_fdr_corrected_list <- list()
# for (DMR in names(shoal_anova_res)) { 
#   #extract the Anova result from the list
#   anova_result <- shoal_anova_res[[DMR]]
#   
#   #extract the uncorrected p-values
#   uncorrected_pvalues <- anova_result$"Pr(>Chisq)"
#   
#   # Perform FDR correction using the Benjamini-Hochberg procedure
#   fdr_corrected_pvalues <- p.adjust(uncorrected_pvalues, method = "hommel")
#   
#   # Store the FDR-corrected p-values in the list
#   shoal_fdr_corrected_list[[DMR]] <- fdr_corrected_pvalues
# }
# 
# #check correlations 
# cor(DMR_all_data$NC_024346.1_10897101, DMR_all_data$NC_024346.1_21149201)
# cor(DMR_all_data$NC_024346.1_10897101, DMR_all_data$NC_024338.1_20758901)
# cor(DMR_all_data$NC_024346.1_21149201, DMR_all_data$NC_024338.1_20758901)
# 
# 
# 
# ## Extra code ##
# ## maybe figure out permutations?? 
# # Number of permutations
# n_permutations <- 100
# 
# # Initialize vectors to store performance metrics
# original_r2 <- c()
# null_r2 <- numeric(n_permutations)
# 
# # Calculate original R^2
# original_r2 <- perf(pls_model)$adjr2
# 
# # Permutation testing
# for (i in 1:n_permutations) {
#   # Permute the response variable
#   permuted_Y <- sample(Y)
#   
#   # Fit null PLSR model with permuted Y
#   null_pls_model <- pls(X, permuted_Y, ncomp = 2)
#   
#   # Calculate R^2 for null model
#   null_r2[i] <- perf(null_pls_model)$adjr2
# }
# 
# # Calculate p-value
# p_value <- mean(null_r2 >= original_r2)
# 
# 
#   
# ## try sPLS2
# #determine the number of dimensions
# tune.pls2.DMR <- pls(X = X, Y = Y, ncomp = 20, mode = "regression")
# Q2.pls2.DMR <- perf(tune.pls2, validation = "Mfold", folds = 10, nrepeat = 5)
# plot(Q2.pls2, criterion = "Q2.total")
# 
# #determine number of variables 
# list.keepX <- c(seq(2, 24, 2))
# list.keepY <- c(1:3)
# 
# tune.spls.all <- tune.spls(X, Y, test.keepX = list.keepX, test.keepY = list.keepY, ncomp = 2,
#                            nrepeat = 1, folds = 10, mode = "regression", measure = "cor")
# 
# plot(tune.spls.all)
# 
# choice.keepX <- tune.spls.all$choice.keepX
# choice.keepY <- tune.spls.all$choice.keepY
# choice.ncomp <- length(choice.keepX)
# 
# #run spls
# spls2.all <- spls(X, Y, ncomp = choice.ncomp, keepX = choice.keepX, keepY = 3, mode = "regression")
# 
# spls2.all$prop_expl_var
# 
# selectVar(spls2.all, comp = 1)$Y
# 
# #check VIP
# vip.all <- vip(spls2.all)
# 
# #look at stability
# perf.spls.all <- perf(spls2.all, validation = "Mfold", folds = 10, nrepeat = 100, progressBar = TRUE)
# 
# stab <- perf.spls.all$features$stability.X
# 
# ## try PCA analysis ##
# library(factoextra)
# 
# pca.DMR <- prcomp(DMR_all_data[1:22])
# 
# pca.DMS <- prcomp(DMS_all_data[1:1500])
# 
# fviz_eig(pca.DMR)
# fviz_eig(pca.DMS)
# 
# res.ind <- get_pca_ind(pca.DMR)
# 
# coord <- res.ind$coord[,1:5]
# 
# pca.behav.dat <- merge(coord, behav.dat, by = "row.names")
# 
# dim.list <- c("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5")
# 
# #run lmms for dist
# pca_dist_lmm_list <- list()
# pca_dist_anova_res <- list()
# 
# for (dim in dim.list) {
#   #run model
#   #run lmm
#   lmm <- lmer(dist_cm ~ get(dim) + sex + cue + (1|tank), REML = TRUE, data = pca.behav.dat)
#   
#   print(summary(lmm))
#   
#   #save model
#   pca_dist_lmm_list[[dim]] <- lmm
#   #run anova
#   anova_res <- Anova(lmm, type = 3)
#   #save anova results
#   pca_dist_anova_res[[dim]] <- anova_res
# }
# 
# #run lmms for shelt
# pca_shelt_lmm_list <- list()
# pca_shelt_anova_res <- list()
# 
# for (dim in dim.list) {
#   #run model
#   #run lmm
#   lmm <- lmer(time_shelter_s ~ get(dim) + sex + cue + (1|tank), REML = TRUE, data = pca.behav.dat)
#   
#   print(summary(lmm))
#   
#   #save model
#   pca_shelt_lmm_list[[dim]] <- lmm
#   #run anova
#   anova_res <- Anova(lmm, type = 3)
#   #save anova results
#   pca_shelt_anova_res[[dim]] <- anova_res
# }
# 
# #run lmms for shoal
# pca_shoal_lmm_list <- list()
# pca_shoal_anova_res <- list()
# 
# for (dim in dim.list) {
#   #run model
#   #run lmm
#   lmm <- lmer(total_diff_s ~ get(dim) + sex + cue + (1|tank), REML = TRUE, data = pca.behav.dat)
#   
#   print(summary(lmm))
#   
#   #save model
#   pca_shoal_lmm_list[[dim]] <- lmm
#   #run anova
#   anova_res <- Anova(lmm, type = 3)
#   #save anova results
#   pca_shoal_anova_res[[dim]] <- anova_res
# }
