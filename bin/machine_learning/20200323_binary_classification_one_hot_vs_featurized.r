# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the one-hot encoded dataset
dat <- read_csv("data/machine_learning/20200322_1095_training_examples_12angstrom_one_hot_encoded_unscaled.csv")

# Set random seed 
set.seed(1234)

# Split into test and training data
dat_split <- initial_split(dat, strata = "is_active")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("id", "is_active")]
x_test <- dat_test[,!colnames(dat_test) %in% c("id", "is_active")]
y_train <- dat_train$is_active
y_test <- dat_test$is_active

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)

# Random forest with one-hot encoding, tuning different mtrys
mtrys <- c(round(log2(ncol(df_train)), 0), round(sqrt(ncol(df_train)), 0), round(ncol(df_train)/2, 0))
mtrys # number of variables available for splitting at each tree node

rf_grid <- expand.grid(mtry = mtrys,
                       splitrule = "gini",
                       min.node.size = 1)
rf <- train(
  x = df_train,
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
rf # Accuracy is 80.5%
#saveRDS(rf, "data/machine_learning/models/20200322_rf_with_one_hot_encoding_unscaled.rds")


rf <- readRDS("data/machine_learning/models/20200322_rf_with_one_hot_encoding_unscaled.rds")

# Train model with best parameters
rf_ml <- ranger(y_train ~., data = form_train, num.trees = 1000, 
                 splitrule = as.character(rf$bestTune$splitrule),
                 mtry = rf$bestTune$mtry, 
                 min.node.size = rf$bestTune$min.node.size,
                 importance = "permutation")

rf_pred <- predict(rf_ml, data = form_test)
cm_rf <- confusionMatrix(rf_pred$predictions, as.factor(dat_test$is_active))
cm_rf

# Now compare to engineered features
dat <- read_csv("data/machine_learning/20200111_1095_training_examples_12angstrom_prot_features.csv")

# Split into test and training data
dat_split <- initial_split(dat, strata = "is_active")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("id", "is_active")]
x_test <- dat_test[,!colnames(dat_test) %in% c("id", "is_active")]
y_train <- dat_train$is_active
y_test <- dat_test$is_active

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)

# Random forest with one-hot encoding
mtrys <- c(round(log2(ncol(df_train)), 0), round(sqrt(ncol(df_train)), 0), round(ncol(df_train)/2, 0))

rf_grid <- expand.grid(mtry = mtrys,
                       splitrule = "gini",
                       min.node.size = 1)
rf2 <- train(
  x = df_train,
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
rf2 # 83.1% training set accuracy

#saveRDS(rf2, "data/machine_learning/models/20200322_rf_with_feature_engineering.rds")
rf2 <- readRDS("data/machine_learning/models/20200322_rf_with_feature_engineering.rds")

# Train model with best parameters
rf2_ml <- ranger(y_train ~., data = form_train, num.trees = 1000, 
                splitrule = as.character(rf2$bestTune$splitrule),
                mtry = rf2$bestTune$mtry, 
                min.node.size = rf2$bestTune$min.node.size,
                importance = "permutation")

rf2_pred <- predict(rf2_ml, data = form_test)
as.factor(dat_test$is_active)
rf2_pred$predictions

cm_rf2 <- confusionMatrix(rf2_pred$predictions, as.factor(dat_test$is_active))
cm_rf2

