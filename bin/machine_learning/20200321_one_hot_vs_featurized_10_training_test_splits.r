# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "readxl", "ggpubr", "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the dataset
dat <- read_csv('data/machine_learning/20200322_1095_training_examples_12angstrom_one_hot_encoded_unscaled.csv')

# Read in the random forest model from cross-validation
rf_20200322 <- readRDS("data/machine_learning/models/20200322_rf_with_one_hot_encoding_unscaled.rds")
rf_20200322$bestTune$mtry
rf_20200322$bestTune$splitrule
rf_20200322$bestTune$min.node.size

# Set random seed 
set.seed(5678)  

rf_models <- vector(mode = "list",
                    length = 10)

# Split randomly into 10 test/training sets
for(i in 1:10) {
  # Split into test and training data
  dat_split <- initial_split(dat, strata = "is_active")
  dat_train <- training(dat_split)
  dat_test  <- testing(dat_split)
  
  # Define variables of interest
  x_train <- dat_train %>%
    dplyr::select(-id, -is_active)
  x_test <- dat_test %>%
    dplyr::select(-id, -is_active)
  
  y_train <- dat_train$is_active
  y_test <- dat_test$is_active
  
  df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)
  
  ## Train model using set parameters (determined by cross-validation)
  rf_1 <- ranger(y_train ~., data = form_train, num.trees = 1000, 
                 splitrule = rf_20200322$bestTune$splitrule,
                 mtry = rf_20200322$bestTune$mtry, 
                 min.node.size = rf_20200322$bestTune$min.node.size,
                 importance = "permutation")
  
  rf_models[[i]]$models <- rf_1
  rf_models[[i]]$oob <- 1 - rf_1$prediction.error # Classifcation accuracy 1 - OOB error
  
  rf_1_pred <- predict(rf_1, form_test)
  
  rf_models[[i]]$confmat <- confusionMatrix(rf_1_pred$predictions, as.factor(dat_test$is_active))
  
  vimp_1 <- data.frame(cbind(sort(rf_1$variable.importance, decreasing = T),
                             names(sort(rf_1$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
    dplyr::rename(importance = X1,
                  variable = X2) %>%
    mutate(importance = as.numeric(importance)) %>%
    dplyr::slice(1:30)
  vimp_1
  rf_models[[i]]$vimp <- vimp_1
  
  vp1 <- ggplot(data = vimp_1,
                aes(x=reorder(variable,importance), y=importance, fill=importance))+
    geom_bar(stat="identity", position="dodge")+ coord_flip()+
    ylab("Variable Importance")+
    xlab("")+
    guides(fill=F)+
    scale_fill_gradient(low="red", high="blue")
  rf_models[[i]]$vimp_plot <- vp1
  
}

mean_training_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$oob }) 
mean_training_vec
mean_testing_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$confmat$overall[1]})
mean_testing_vec

dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) 
sd_df <- dat_df %>%
  dplyr::summarise_each(sd)
sd_df
colMeans(dat_df)
