# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Set random seed 
set.seed(1234)

# Read in the full dataset
dat <- read_csv("data/machine_learning/20200111_1095_training_examples_12angstrom_prot_features.csv")
colnames(dat)
# Data 
x_train <- dat[,!colnames(dat) %in% c("id", "is_active")]
y_train <- dat$is_active
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat$id)

# Random forest with one-hot encoding, tuning different mtrys
mtrys <- c(round(log2(ncol(dat)), 0), round(sqrt(ncol(dat)), 0), round(ncol(dat)/2, 0))
mtrys # number of variables available for splitting at each tree node

rf_grid <- expand.grid(mtry = mtrys,
                       splitrule = "gini",
                       min.node.size = 1)
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
rf # Accuracy is 80.6%

rf$results # 82.4% accuracy

# Now train on the entire dataset for web app model
rf_full_ss_noxval <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                            mtry = round(ncol(dat)/2, 0), min.node.size = 1,
                            importance = "permutation", probability = TRUE)
rf_full_ss_noxval # OOB error is 0.1218886, or about 87.8% prediction accruacy


saveRDS(rf_full_ss_noxval, "data/20200323_rf_fullset_ss_noxval.rds")
 
#saveRDS(rf, "data/machine_learning/models/20200322_rf_with_one_hot_encoding_unscaled.rds")

   