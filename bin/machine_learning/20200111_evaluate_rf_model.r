# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "gbm","nnet","e1071","svmpath","lars",
               "glmnet", "svmpath", "readxl", "ggpubr", "doMC", "doParallel")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the random forest model
rf_20200111 <- readRDS("data/machine_learning/models/20200111_rf_10foldcv.rds")
rf_20200111$bestTune
getTrainPerf(rf_20200111)
