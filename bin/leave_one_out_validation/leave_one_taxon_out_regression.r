# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", "readxl", "ggpubr",
               "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the data
dat <- read_csv("data/machine_learning/20200111_1095_training_examples_12angstrom_prot_features_activity_all.csv") %>%
  dplyr::mutate(spec = word(id, sep = "_", 1)) %>%
  dplyr::filter(activity > 0)
dat$spec

# Read in the taxonomic identification
rawtax <- read_excel("data/OleA_taxonomic_classification.xlsx") 
xantho <- rawtax[grep("Xanthomonas", rawtax$organism),]
xantho$organism <- c("Xanthomonas_campestris_OleA")
xantho$genus <- "Xanthomonas_campestris"

tax <- rawtax  %>%
  bind_rows(xantho)

tax$organism <- gsub("_", " ", tax$organism)
dtax <- tax %>%
  dplyr::mutate(spec = paste0(word(organism, sep = " ", 1), " ", word(organism, sep = " ", 2)))

# Fix discrepancies in merging names
pseudo1 <- dtax$spec[grep("Pseudoxanthomonas", dtax$spec)]
pseudo2 <- dat$spec[grep("Pseudoxanthomonas", dat$spec)] 
dat$spec[grep("Pseudoxanthomonas", dat$spec)] <- pseudo1

leif1 <- dtax$spec[grep("Leifsonia", dtax$spec)]
leif2 <- dat$spec[grep("Leifsonia", dat$spec)] 
dat$spec[grep("Leifsonia", dat$spec)] <- leif1

# Merge with the tax df
merg <- dtax %>%
  inner_join(dat, by = "spec") %>%
  dplyr::select(-organism, -genus, -family, -order, -phylum, -spec)
uniq_tax <- unique(merg$class)

# Set random seed 
set.seed(5678)  

# Create a list of models
rf_models <- vector(mode = "list",
                    length = length(uniq_tax))

# Leave-one-taxon-out validation
for(i in 1:length(uniq_tax)) {
  
  # Identify compound to be 'left out' of model training
  loo_tax <- uniq_tax[i]
  
  # Leave one taxon out
  dat_train <- merg[merg$class != loo_tax,]
  dat_test  <- merg[merg$class == loo_tax,]
  
  # Define variables of interest
  x_train <- dat_train %>%
    dplyr::select(-id, -activity, -class)
  x_test <- dat_test %>%
    dplyr::select(-id, -activity, -class)
  
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

accuracy <- round(mean_testing_vec * 100, 2)
rf_models[[which.min(mean_testing_vec)]] 

test_set_size <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$test_set_size})
train_set_size <- 500 - test_set_size

dat_df <- tibble(cmpnds[1:14]) %>%
  bind_cols(tibble(accuracy)) %>%
  bind_cols(tibble(test_set_size)) %>%
  bind_cols(tibble(train_set_size)) %>%
  arrange(accuracy)
dat_df

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

dat_df <- tibble(uniq_tax) %>%
  bind_cols(tibble(accuracy)) %>%
  bind_cols(tibble(test_set_size)) %>%
  bind_cols(tibble(train_set_size)) %>%
  arrange(accuracy)
dat_df

colnames(dat_df) <- c("left-out taxon", "leave-one-taxon-out RMSE")
write_csv(dat_df, "output/Leave-one-taxon-out-regression-analysis.csv")





