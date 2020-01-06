# Read in the raw data
# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest","gbm","nnet","e1071","svmpath","lars",
               "glmnet","svmpath", "readxl")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the molecular features
molec_fts <- read_csv("data/machine_learning/PC7_molecular_descriptors.csv") %>%
  dplyr::rename(substrate = V1)
molec_fts$substrate <- gsub(" ", "\\.", molec_fts$substrate)
molec_fts$substrate[molec_fts$substrate == "oxidazole"] <- "oxadiazole"
head(molec_fts)

# Read in the activity data
activity <- read_csv("data/machine_learning/20191218_all_cmpnds_avg_log_slopes_for_modeling.csv")

# Read in the cavity volumes
cavity_fts <- read_excel("data/caver_models/CastP_volume_rawdata.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::mutate(genus = gsub("\\.pdb", "", word(filename, sep = "_", 4))) %>%
  dplyr::mutate(species = gsub("\\.pdb", "", word(filename, sep = "_", 5))) %>%
  dplyr::mutate(org = paste0(genus, " ", species)) 

# Now merge everything...
comb <- activity %>%
  dplyr::left_join(., molec_fts, by = "substrate") %>%
  dplyr::left_join(., cavity_fts, by = "org") %>%
  dplyr::mutate(id = paste0(org, "_", substrate)) %>%
  dplyr::mutate(is_active = as.factor(case_when(activity == 0 ~ "N",
                                                TRUE ~ "Y"))) %>%
  dplyr::select(id, contains("sa"),  contains("PC"), is_active)

# Now check for missing data and duplicates
dedup <- comb[complete.cases(comb),] 
dedup <- dedup[!duplicated(comb),]

# Only keep variables with nonzero variance
nozdat <- nearZeroVar(dedup, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE]

# Set random seed 
set.seed(20201002)

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
table(dat$is_active)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)

# Quick  test

rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  # tuneGrid = tgrid,
  num.trees = 500,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  preProcess = c("center", "scale"),
  importance = "permutation") 

# Confusion matrix
getTrainPerf(rf) # Training set accuracy is 82% for random forest

# Try prediction
rf$bestTune$splitrule
rf$bestTune$mtry
ncol(form_train)

rf_ml <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini", #splitrule = as.character(rf$bestTune$splitrule),
                mtry = sqrt(ncol(form_train)), #rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
                importance = "permutation")

rf_pred <- predict(rf, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$is_active))
cm_rf # 80.22% test set accuracy

sink("output/machine_learning/random_forest_12angstroms_binary_classification_results_20200201.txt")
cm_rf
sink()

vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T),
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:15)
vimp

