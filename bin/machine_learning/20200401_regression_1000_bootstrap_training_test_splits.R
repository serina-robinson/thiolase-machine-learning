# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "readxl", "ggpubr", "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the dataset
dat <- read_csv('data/machine_learning/20200111_550_regression_training_examples_12angstrom_prot_features.csv')

# Read in the random forest model from cross-validation
rf_20200322 <- readRDS("data/machine_learning/models/20200104_final_rf_regression_model.rds")

# mtry splitrule min.node.size
# 168  variance             5
rf_20200322$mtry # 168
rf_20200322$splitrule # variance
rf_20200322$min.node.size #5

# Set random seed 
set.seed(5678)  

rf_models <- vector(mode = "list",
                    length = 1000)

# Split randomly into 10 test/training sets
for(i in 1:1000) {
  # Split into test and training data
  dat_split <- initial_split(dat)
  dat_train <- training(dat_split)
  dat_test  <- testing(dat_split)
  
  # Define variables of interest
  x_train <- dat_train %>%
    dplyr::select(-id, -activity)
  x_test <- dat_test %>%
    dplyr::select(-id, -activity)
  
  y_train <- dat_train$activity
  y_test <- dat_test$activity
  
  df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)
  
  ## Train model using set parameters (determined by cross-validation)
  rf_1 <- ranger(y_train ~., data = form_train, 
                 num.trees = 1000, 
                 splitrule = "variance",
                 mtry = 168, 
                 min.node.size = 5,
                 importance = "permutation")
  
  rf_models[[i]]$models <- rf_1
  rf_models[[i]]$train_rmse <- sqrt(rf_1$prediction.error) # Classifcation accuracy 1 - OOB error
  
  rf_1_pred <- predict(rf_1, form_test)
  
  rf_models[[i]]$test_rmse <- Metrics::rmse(rf_1_pred$predictions, y_test)
}

mean_training_vec <- unlist(sapply(1:length(rf_models), function(x) { rf_models[[x]]$train_rmse })) 
traindf <- tibble(mean_training_vec)
write_tsv(traindf,"output/regression_mean_training_accuracy_1000_test_train_splits.txt")

mean_testing_vec <- unlist(sapply(1:length(rf_models), function(x) { rf_models[[x]]$test_rmse}))
testdf <- tibble(mean_testing_vec)
write_tsv(testdf,"output/regression_mean_testing_accuracy_1000_test_train_splits.txt")

p1 <- ggplot(traindf) +
  geom_histogram(aes(x = mean_training_vec, y = ..density..), alpha = 0.7,
                 binwidth = 0.005, fill = "gray80", color = "black") +
  theme_pubr(base_size = 17) +
  xlab("Mean training set accuracy") +
  ylab("Frequency")
p1


p2 <- ggplot(testdf) +
  geom_histogram(aes( x = mean_testing_vec, y = ..density..), alpha = 0.7,
                 binwidth = 0.005, fill = "gray80", color = "black") +
  theme_pubr(base_size = 17) +
  xlab("Mean testing set accuracy") +
  ylab("Frequency")
p2


# Combine into a single histogram
traindf$split <- "training"
colnames(traindf) <- c("score", "split")
trainsumm <- data.frame(tibble(t(summary(traindf$score))), stringsAsFactors = F)
trainsumm

testdf$split <- "testing"
colnames(testdf) <- c("score", "split")
summary(testdf$score)
testsumm <- data.frame(tibble(t(summary(testdf$score))), stringsAsFactors = F)
testsumm


allsumm <- data.frame(c(trainsumm,testsumm)) 
write_csv(allsumm, "output/regression_train_test_1000_summaries.csv")


alldf <- traindf %>%
  bind_rows(testdf)
head(alldf)
alldf$split <- as.factor(alldf$split)

pdf("output/1000_regression_training_testing_split_hist.pdf", width = 7, height = 4)
p3 <- ggplot(alldf, aes(x= score, fill = split)) +
  geom_histogram(alpha = 0.5,
                 binwidth = 0.005, position = 'identity') +
  theme_pubr(base_size = 12) +
  xlab("Root-mean square error") +
  ylab("Count") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
p3
dev.off()

dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) 
colMeans(dat_df)

sd_df <- dat_df %>%
  dplyr::summarise_each(sd)
sd_df



