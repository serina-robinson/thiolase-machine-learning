# Read in the raw data
# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "gbm","nnet","e1071","svmpath","lars",
               "glmnet", "svmpath", "readxl")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the principal componenets of molecular features
molec_fts <- read_csv("data/machine_learning/PC7_molecular_descriptors.csv") %>%
  dplyr::rename(substrate = V1)
molec_fts$substrate <- gsub(" ", "\\.", molec_fts$substrate)
molec_fts$substrate[molec_fts$substrate == "oxidazole"] <- "oxadiazole"

# Read in the raw molecular features
chem_descr <- read_csv("data/substrate_comparisons/15pNPs_159_selected_molecular_properties.csv") %>%
  #dplyr::select(cmpnd_abbrev, nB, MW, AROMATIC, O, N, Cl, MLogP, nRotB, VABC, nAtomLAC) %>%
  dplyr::mutate(substrate = gsub(" ", "\\.", cmpnd_abbrev)) %>%
  dplyr::select(-cmpnd_abbrev)
chem_descr$substrate[chem_descr$substrate == "oxidazole"] <- "oxadiazole"

# Read in the sequence features 
seq_fts <- read_csv("data/machine_learning/73_12angstrom_4aa_features.csv") %>%
  dplyr::mutate(raw_org = word(enzyme, sep = "\\.1", 2)) %>%
  dplyr::mutate(org = paste0(word(raw_org, sep = "_", 2), " ", word(raw_org, sep = "_", 3))) %>%
  dplyr::select(-raw_org) # remove enzyme

seq_fts$org[seq_fts$enzyme == "4KU5_Xanthomonas_campestris"] <- "Xanthomonas campestris"

# Read in the activity data
activity <- read_csv("data/machine_learning/20191218_all_cmpnds_avg_log_slopes_for_modeling.csv")

# Fix discrepancies in merging names
pseudo1 <- seq_fts$org[grep("Pseudoxanthomonas", seq_fts$org)]
pseudo2 <- activity$org[grep("Pseudoxanthomonas", activity$org)][1]
seq_fts$org[grep("Pseudoxanthomonas", seq_fts$org)] <- pseudo2

leif1 <- seq_fts$org[grep("Leifsonia", seq_fts$org)]
leif2 <- activity$org[grep("Leifsonia", activity$org)][1]
seq_fts$org[grep("Leifsonia", seq_fts$org)] <- leif2

# Now merge everything...
comb <- activity %>%
  dplyr::left_join(., molec_fts, by = "substrate") %>%
  dplyr::left_join(., chem_descr, by = "substrate") %>%
  #dplyr::left_join(., cavity_fts, by = "org") %>%
  dplyr::left_join(., seq_fts, by = "org")

# Now remove duplicate rows (hopefully there aren't any)
dedup <- comb[complete.cases(comb),] # no duplicates
dedup <- dedup[!duplicated(dedup),]

# Only keep variables with nonzero variance
nozdat <- nearZeroVar(dedup, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE]
which_rem

# write_csv(dedup, "data/machine_learning/20191228_1095_training_examples_12angstrom_features.csv")
dat <- dedup %>%
  dplyr::mutate(id = paste0(org, "_", substrate)) %>%
  dplyr::mutate(is_active = as.factor(case_when(activity == 0 ~ "N",
                                                TRUE ~ "Y"))) %>%
  dplyr::select(-which_rem, -org, -enzyme, -activity, -substrate, -IUPAC, -SMILES, -substrate)
colnames(dat)

# Set random seed 
set.seed(1234)

# Split into test and training data
dat_split <- initial_split(dat, strata = "is_active")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %

##### Test RF with only raw chemical descriptors and sequence features

# Define variables of interest
x_train <- dat_train %>%
  dplyr::select(-id, -is_active)
colnames(x_train)

x_test <- dat_test %>%
  dplyr::select(-id, -is_active)
colnames(x_train)
colnames(x_test)

y_train <- dat_train$is_active
y_test <- dat_test$is_active

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)

# rf_grid <- expand.grid(mtry = c(25, 50, 75, 100, 125, 150, 175, 200),
#                        splitrule = "extratrees",
#                        min.node.size = 1)

# Quick  test
rf_1 <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 500,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  preProcess = c("center", "scale"),
  importance = "permutation") 

saveRDS(rf_1, "data/machine_learning/models/20200106_rf_all_chem_descr.rds")

# Confusion matrix
getTrainPerf(rf_1) # 80.6% test set accuracy

# Try prediction
rf_1$bestTune$splitrule
rf_1$bestTune$mtry
rf_1$bestTune$min.node.size



rf_ml <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = as.character(rf_1$bestTune$splitrule),
                mtry = rf_1$bestTune$mtry, min.node.size = rf_1$bestTune$min.node.size,
                importance = "permutation", probability = T)

rf_pred <- predict(rf_1, newdata = form_test)

cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$is_active))
cm_rf # 80.22% test set accuracy


sink("output/machine_learning/20200106_rf_all_descriptors_classification_results.txt")
cm_rf
sink()

vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T),
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:100)
vimp
