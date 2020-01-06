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
  dplyr::select(cmpnd_abbrev, nB, MW, AROMATIC, O, N, Cl, MLogP, nRotB, VABC, nAtomLAC) %>%
  dplyr::mutate(substrate = gsub(" ", "\\.", cmpnd_abbrev)) %>%
  dplyr::select(-cmpnd_abbrev)
chem_descr$substrate[chem_descr$substrate == "oxidazole"] <- "oxadiazole"

# Read in the cavity volumes
cavity_fts <- read_excel("data/caver_models/CastP_volume_rawdata.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::mutate(genus = gsub("\\.pdb", "", word(filename, sep = "_", 4))) %>%
  dplyr::mutate(species = gsub("\\.pdb", "", word(filename, sep = "_", 5))) %>%
  dplyr::mutate(org = paste0(genus, " ", species)) %>%
  dplyr::select(org, volume_sa, area_sa)

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
  dplyr::left_join(., cavity_fts, by = "org") %>%
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
  dplyr::select(-which_rem, -org, -substrate) %>%
  dplyr::filter(is_active == "Y")

hist(dat$activity)

# Set random seed 
set.seed(20200105)

# Split into test and training data
dat_split <- initial_split(dat)
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %

##### Test RF with only raw chemical descriptors and sequence features

# Define variables of interest
x_train <- dat_train %>%
  dplyr::select(-contains("PC"), -volume_sa, -area_sa, -enzyme, -activity, -id, -is_active)
x_test <- dat_test %>%
  dplyr::select(-contains("PC"), -volume_sa, -area_sa, -enzyme, -activity, -id, -is_active)
colnames(x_train)
colnames(x_test)

y_train <- dat_train$activity
y_test <- dat_test$activity

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)



rf_grid <- expand.grid(mtry = c(25, 50, 75, 100, 125, 150, 175, 200),
                       splitrule = "extratrees",
                       min.node.size = 1)

# Quick  test
rf_1 <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 5,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  tuneGrid = rf_grid,
  num.trees = 1000,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  preProcess = c("center", "scale"),
  importance = "permutation") 

# Confusion matrix
getTrainPerf(rf_1) # RMSE is 0.257
rmse
# Try prediction
rf_ml_1 <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = as.character(rf_1$bestTune$splitrule),
                   mtry = rf_1$bestTune$mtry, min.node.size = rf_1$bestTune$min.node.size,
                   importance = "permutation")

rf_1_pred <- predict(rf_1, newdata = form_test)

rmse(rf_1_pred, y_test)
Metrics::rmse(y_test, rf_1_pred) # 0.274


#cm_rf_1 <- confusionMatrix(rf_1_pred, as.factor(dat_test$is_active))
#cm_rf_1 # 76.6% test classification accuracy with limited chemical features

# sink("output/machine_learning/20200104_rf1_rawchem_and_seq_regression_results.txt")
# cm_rf_1
# sink()

vimp_1 <- data.frame(cbind(sort(rf_ml_1$variable.importance, decreasing = T),
                         names(sort(rf_ml_1$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:30)
vimp_1

pdf("output/machine_learning/rf1_chem_descr_only_regression_vimp.pdf", width = 6, height = 6)
vp1 <- ggplot(data = vimp_1,
       aes(x=reorder(variable,importance), y=importance, fill=importance))+
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
vp1
dev.off()

# Try training a mars model with the chemical descriptors
mars_grid <- expand.grid(degree = 1:2, nprune = seq(2, 26, by = 2))

mars_mod <- train(
  x = df_train,
  y = y_train,
  method = "earth",
  tuneGrid = mars_grid,
  preProcess = c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", 
                           repeats = 5, 
                           number = 10,
                           savePredictions = "final"))

# Training and test performance
getTrainPerf(mars_mod) # RMSE is 0.274
mars_mod$results
mars_mod$finalModel
mars_pred <- predict(mars_mod, newdata = form_test)
mars_pred

mars_mod$finalModel$coefficients
mars_imp <- varImp(mars_mod)
ggplot(mars_imp, top = 20) + xlab("")
Metrics::rmse(y_test, mars_pred) # 0.32

##### RF with only chem PCs and sequence feature

# Define variables of interest
x_train <- dat_train %>%
  dplyr::select(contains("PC"), contains("4KU5"))
x_test <- dat_test %>%
  dplyr::select(contains("PC"), contains("4KU5"))

y_train <- dat_train$activity
y_test <- dat_test$activity

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)

# Quick  test
rf_2 <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 5,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  tuneGrid = rf_grid,
  num.trees = 1000,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  preProcess = c("center", "scale"),
  importance = "permutation") 

# Confusion matrix
getTrainPerf(rf_2) # Training set accuracy is 82.9% for random forest

# Try prediction
rf_2$results
rf_2$bestTune$splitrule
rf_2$bestTune$mtry # mtry = 125
rf_2$bestTune$min.node.size

rf_ml_2 <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = as.character(rf_2$bestTune$splitrule),
                  mtry = rf_2$bestTune$mtry, min.node.size = rf_2$bestTune$min.node.size,
                  importance = "permutation")

rf_2_pred <- predict(rf_2, newdata = form_test)
Metrics::rmse(y_test, rf_2_pred) # 0.258773 RMSE

# cm_rf_2 <- confusionMatrix(rf_2_pred, as.factor(dat_test$is_active))
# cm_rf_2 # 79.5% test accuracy


channel_a <- c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353)
channel_b <- c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246)
channels <- c(channel_a, channel_b)

rf_ml_2$variable.importance != 0
imps <- word(names(rf_ml_2$variable.importance)[rf_ml_2$variable.importance != 0], -1, sep = "_")
imps
imps <- as.numeric(imps[!grepl("PC", imps)])
imps[imps %in% channel_a]
imps[imps %in% channel_b]


# sink("output/machine_learning/20200104_rf2_PC_and_seq_regression_classification_results.txt")
# cm_rf_2
# sink()

vimp_2 <- data.frame(cbind(sort(rf_ml_2$variable.importance, decreasing = T),
                           names(sort(rf_ml_2$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:30)
vimp_2

pdf("output/machine_learning/rf2_PC_seqfeats_regression_classification_vimp.pdf", width = 6, height = 6)
vp2 <- ggplot(data = vimp_2,
       aes(x=reorder(variable,importance), y=importance, fill=importance))+
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
vp2
dev.off()

##### RF with chem PCs, sequence features, and surface areas

# Define variables of interest
x_train <- dat_train %>%
  dplyr::select(contains("PC"), contains("4KU5"), contains("_sa"))
x_test <- dat_test %>%
  dplyr::select(contains("PC"), contains("4KU5"), contains("_sa"))

y_train <- dat_train$activity
y_test <- dat_test$activity

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)

# Quick  test
rf_3 <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 5,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  tuneGrid = rf_grid,
  num.trees = 1000,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  preProcess = c("center", "scale"),
  importance = "permutation") 

# Confusion matrix
getTrainPerf(rf_3) # Training set RMSE is 0.244
rf_3
rf_3$bestTune$splitrule
rf_3$bestTune$mtry
rf_3$bestTune$min.node.size

# Try prediction
rf_ml_3 <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = as.character(rf_3$bestTune$splitrule),
                  mtry = rf_3$bestTune$mtry, min.node.size = rf_3$bestTune$min.node.size,
                  importance = "permutation")

rf_3_pred <- predict(rf_3, newdata = form_test)
Metrics::rmse(rf_3_pred, y_test) # 0.258 RMSE testing

# cm_rf_3 <- confusionMatrix(rf_3_pred, as.factor(dat_test$is_active))
# cm_rf_3 # 79.5% test accuracy

# sink("output/machine_learning/20200104_rf3_PC_and_seq_and_volume_regression_results.txt")
# cm_rf_3
# sink()

vimp_3 <- data.frame(cbind(sort(rf_ml_3$variable.importance, decreasing = T),
                           names(sort(rf_ml_3$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:30)
vimp_3

pdf("output/machine_learning/rf3_PC_seqfeats_sa_regression_vimp.pdf", width = 6, height = 6)
vp3 <- ggplot(data = vimp_3,
       aes(x=reorder(variable,importance), y=importance, fill=importance))+
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
vp3 
dev.off()

pdf("output/machine_learning/20200104_regression_comparisons.pdf")
plot_grid(vp1, vp2, vp3, ggplot(mars_imp, top = 20), ncol = 4)
dev.off()






