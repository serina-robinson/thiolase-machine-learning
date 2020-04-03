# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "readxl", "ggpubr", "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the dataset
dat <- read_csv('data/machine_learning/20200111_1095_training_examples_12angstrom_prot_features.csv')

# Read in the random forest model from cross-validation
rf_20200322 <- readRDS("data/machine_learning/models/20200111_rf_10foldcv.rds")
rf_20200322$bestTune$mtry
rf_20200322$bestTune$splitrule
rf_20200322$bestTune$min.node.size

# Set random seed 
set.seed(5678)  

rf_models <- vector(mode = "list",
                    length = 1000)

# Split randomly into 10 test/training sets
for(i in 1:1000) {
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
  rf_1 <- ranger(y_train ~., data = form_train, 
                 num.trees = 1000, 
                 splitrule = "extratrees",
                 mtry = 334, 
                 min.node.size = 1,
                 importance = "permutation")
  
  rf_models[[i]]$models <- rf_1
  rf_models[[i]]$oob <- 1 - rf_1$prediction.error # Classifcation accuracy 1 - OOB error
  
  rf_1_pred <- predict(rf_1, form_test)
  
  rf_models[[i]]$confmat <- confusionMatrix(rf_1_pred$predictions, as.factor(dat_test$is_active))
}

mean_training_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$oob }) 
traindf <- tibble(mean_training_vec)
write_tsv(traindf,"output/classification_mean_training_accuracy_1000_test_train_splits.txt")

p1 <- ggplot(traindf) +
  geom_histogram(aes( x = mean_training_vec, y = ..density..), alpha = 0.7,
                 binwidth = 0.005, fill = "gray80", color = "black") +
  theme_pubr(base_size = 17) +
  xlab("Mean training set accuracy") +
  ylab("Frequency")
p1

mean_testing_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$confmat$overall[1]})
testdf <- tibble(mean_testing_vec)

p2 <- ggplot(testdf) +
  geom_histogram(aes( x = mean_testing_vec, y = ..density..), alpha = 0.7,
                 binwidth = 0.005, fill = "gray80", color = "black") +
  theme_pubr(base_size = 17) +
  xlab("Mean testing set accuracy") +
  ylab("Frequency")
p2

write_tsv(tibble(mean_testing_vec),"output/classification_mean_testing_accuracy_1000_test_train_splits.txt")

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
write_csv(allsumm, "output/classification_train_test_1000_summaries.csv")


alldf <- traindf %>%
  bind_rows(testdf)
head(alldf)
alldf$split <- as.factor(alldf$split)

pdf("output/1000_binary_classification_training_testing_split_hist.pdf", width = 7, height = 4)
p3 <- ggplot(alldf, aes(x= score, fill = split)) +
  geom_histogram(alpha = 0.5,
                 binwidth = 0.005, position = 'identity') +
  theme_pubr(base_size = 12) +
  xlab("Classification accuracy") +
  ylab("Count") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
p3
dev.off()


dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) 
dat_df
colMeans(dat_df)

sd_df <- dat_df %>%
  dplyr::summarise_each(sd)
sd_df



