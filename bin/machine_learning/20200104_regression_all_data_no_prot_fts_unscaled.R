# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Set random seed 
set.seed(1234)

# Read in the full dataset
dat <- read_csv("data/machine_learning/20200323_550_regression_training_examples_no_prot_feats_unscaled.csv")

# Data 
x_train <- dat[,!colnames(dat) %in% c("id", "activity")]
y_train <- dat$activity
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat$id)
y_train

# Random forest with one-hot encoding, tuning different mtrys
mtrys <- c(round(log2(ncol(dat)), 0), round(sqrt(ncol(dat)), 0), round(ncol(dat)/2, 0), 
           round(ncol(dat) * (2/3), 0), round(ncol(dat) * (5/6), 0))
mtrys # number of variables available for splitting at each tree node

rf_grid <- expand.grid(mtry = mtrys,
                       splitrule = c("variance", "extratrees"),
                       min.node.size = c(1, 3, 5))
rf <- train(
  x = x_train,
  y = y_train,
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final",
                           returnResamp = "all"),
  num.trees = 1000,
  verbose = TRUE,
  importance = "permutation") 


getTrainPerf(rf) # RMSE 0.242
rf$bestTune
# mtry = 140, splitrule = extratrees, min.node.size = 5
saveRDS(rf, "data/machine_learning/models/20200401_rf_regression_all_data_unscaled_xval.rds")

# Now train a simple model using best tuning parameters from 10fold_xval
rf_bin_class_all <- ranger(y_train ~., data = form_train, num.trees = 1000, 
                           splitrule = rf$bestTune$splitrule,
                           mtry = rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
                           importance = "permutation")
rf_bin_class_all 
saveRDS(rf_bin_class_all, "shiny_app/data/models/20200401_rf_regression_all_data_no_prot_unscaled_optimized.rds")
# head(dat)
#saveRDS(rf, "data/machine_learning/models/20200322_rf_with_one_hot_encoding_unscaled.rds")

