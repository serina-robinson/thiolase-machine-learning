# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", "readxl", "ggpubr",
               "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the data
dat <- read_csv("data/machine_learning/20200111_1095_training_examples_12angstrom_prot_features_activity_all.csv") %>%
  mutate(cmpnd = (word(id, sep = "_", -1))) %>%
  dplyr::filter(activity > 0)

cmpnds <- unique(word(dat$id, sep = "_", -1)) # 14 compounds

# Set random seed 
set.seed(5678)  

# Create a list of models
rf_models <- vector(mode = "list",
                    length = length(cmpnds))

# Leave-one-compound-out validation
for(i in 1:length(cmpnds)) {
  
  # Identify compound to be 'left out' of model training
  loo_cmpnd <- cmpnds[i]
  
  # Leave one compound out
  dat_train <- dat[dat$cmpnd != loo_cmpnd,]
  dat_test  <- dat[dat$cmpnd == loo_cmpnd,]
  print(dim(dat_test))
  
  # Define variables of interest
  x_train <- dat_train %>%
    dplyr::select(-id, -activity, -cmpnd)
  x_test <- dat_test %>%
    dplyr::select(-id, -activity, -cmpnd)
  
  y_train <- dat_train$activity
  y_test <- dat_test$activity
  
  df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)
  
  ## Train model using default parameters (gini split-rule, sqrt(m) for mtry)
  rf_1 <- ranger(y_train ~., data = form_train, num.trees = 1000, 
                 mtry = sqrt(ncol(form_train)),
                 importance = "permutation")
  
  rf_models[[i]]$models <- rf_1
  rf_models[[i]]$train_rmse <- sqrt(rf_1$prediction.error) # Classifcation accuracy 1 - OOB error
  
  rf_1_pred <- predict(rf_1, form_test)
  rf_models[[i]]$test_set_size <- nrow(dat_test)
  rf_models[[i]]$test_rmse <- Metrics::rmse(rf_1_pred$predictions, y_test)
}

# Benzoate causes an error because the model has not been trained on that type of chemical features previously

# Will report the results for the 14 other compounds

mean_training_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$train_rmse }) 
mean_training_vec

mean_testing_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$test_rmse})
mean_testing_vec

accuracy <- round(mean_testing_vec, 2)
rf_models[[which.min(mean_testing_vec)]] 

test_set_size <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$test_set_size})
train_set_size <- 500 - test_set_size

dat_df <- tibble(cmpnds[1:14]) %>%
  bind_cols(tibble(accuracy)) %>%
  bind_cols(tibble(test_set_size)) %>%
  bind_cols(tibble(train_set_size)) %>%
  arrange(accuracy)
dat_df

colnames(dat_df) <- c("left-out compound", "leave-one-compound-out RMSE", "Testing set size", "Training set size")
write_csv(dat_df, "output/Leave-one-compound-out-regression-analysis.csv")





