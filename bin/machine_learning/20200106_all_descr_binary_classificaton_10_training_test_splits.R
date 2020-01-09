# Read in the raw data
# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "gbm","nnet","e1071","svmpath","lars",
               "glmnet", "svmpath", "readxl", "ggpubr")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the principal componenets of molecular features
# molec_fts <- read_csv("data/machine_learning/PC7_molecular_descriptors.csv") %>%
#   dplyr::rename(substrate = V1)
# molec_fts$substrate <- gsub(" ", "\\.", molec_fts$substrate)
# molec_fts$substrate[molec_fts$substrate == "oxidazole"] <- "oxadiazole"

# Read in the chemical features
chem_descr <- read_csv("data/substrate_comparisons/15pNPs_159_selected_molecular_properties.csv") %>%
  dplyr::mutate(substrate = gsub(" ", "\\.", cmpnd_abbrev))
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
  dplyr::left_join(., chem_descr, by = "substrate") %>%
  dplyr::left_join(., seq_fts, by = "org")
  
# Now remove duplicate rows (hopefully there aren't any)
dedup <- comb[complete.cases(comb),] # no duplicates
dedup <- dedup[!duplicated(dedup),]

# Only keep variables with nonzero variance
nozdat <- nearZeroVar(dedup, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE]


# write_csv(dedup, "data/machine_learning/20191228_1095_training_examples_12angstrom_features.csv")
dat <- dedup %>%
  dplyr::mutate(id = paste0(org, "_", substrate)) %>%
  dplyr::mutate(is_active = as.factor(case_when(activity == 0 ~ "N",
                                                TRUE ~ "Y"))) %>%
  dplyr::select(-which_rem, -org, -substrate, -enzyme, -activity, -IUPAC, -SMILES, -cmpnd_abbrev) 
colnames(dat)

# Set random seed 
set.seed(1234)

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
                 splitrule = "gini",
                 mtry = ncol(x_train)/2, min.node.size = 1,
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

rf_models[[which.max(mean_testing_vec)]] # 82.42 % testing accuracy
# saveRDS(rf_models[[9]]$models, "data/machine_learning/models/final_random_forest_classification_model.rds")


dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) 
sd_df <- dat_df %>%
  dplyr::summarise_each(sd)
sd_df
colMeans(dat_df)


rf_test_pred <- predict(rf_models[[9]]$models, form_test)

comp_df <- data.frame(cbind(dat_test$id, rf_test_pred$predictions, y_test)) %>%
  dplyr::mutate(y_pred = as.numeric(V2)) %>%
  dplyr::mutate(y_test = as.numeric(y_test)) %>%
  dplyr::mutate(incorrect = abs(y_pred - y_test)) %>%
  dplyr::select(V1, y_test, y_pred, incorrect) %>%
  dplyr::filter(incorrect == 1)
comp_df

rf_models[[1]]$vimp[,1]

### VIMPS
vimp_combined <- lapply(1:length(rf_models), function(x) { rf_models[[x]]$vimp }) %>%
  map_dfr(~ .) %>%
  group_by(variable) %>%
  summarise_at(vars(importance), mean) %>%
  arrange(desc(importance)) %>%
  dplyr::mutate(var_fix = case_when(grepl("PC", variable) ~ variable,
                TRUE ~ paste0(word(variable, -1, sep = "_"), " ", word(variable, 1, sep = "_"))))
vimp_combined


vpavg <- ggplot(data = vimp_combined[1:30,],
              aes(x=reorder(var_fix, importance), y=importance, fill=importance))+
  geom_bar(stat="identity", position="dodge") + coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F) +
  theme_pubr()
  
pdf("output/machine_learning/20200104_rf_classification_10iterations_varimp_avged.pdf", height = 6, width = 5)
vpavg
dev.off()

dev.off()

channel_a <- c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353)
channel_b <- c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246)
channels <- c(channel_a, channel_b)

imps <- stringr::word(vimp_combined$variable, -1, sep = "_")
imps
paste0(unique(imps), collapse = "+")

imps <- as.numeric(imps[!grepl("PC", imps)])
sort(table(imps), decreasing = T)

chana <- imps[imps %in% channel_a]
chanb <- imps[imps %in% channel_b]

table(chana)
table(chanb)


vimp_for_pymol <- lapply(1:length(rf_models), function(x) { rf_models[[x]]$vimp }) %>%
  map_dfr(~ .) %>%
  group_by(variable) %>%
  summarise_at(vars(importance), mean) %>%
  arrange(desc(importance)) %>%
  dplyr::slice(1:30) %>%
  dplyr::filter(!grepl("PC", variable)) 

vimp_for_pymol
imps <- stringr::word(vimp_for_pymol$variable, -1, sep = "_")
chana <- imps[imps %in% channel_a]
chanb <- imps[imps %in% channel_b]
chans <- c(chana, chanb)

no_chans <- imps[!imps %in% chans]
no_chans
paste0(unique(no_chans), collapse = "+")
paste0(unique(chans), collapse = "+")
#### Need to make an ROC curve






# channel_a <- c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353)
# channel_b <- c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246)
# channels <- c(channel_a, channel_b)
# 
# rf_ml_2$variable.importance != 0
# imps <- word(names(rf_ml_2$variable.importance)[rf_ml_2$variable.importance != 0], -1, sep = "_")
# imps <- as.numeric(imps[!grepl("PC", imps)])
# imps[imps %in% channel_a]
# imps[imps %in% channel_b]
# 
# 
# 
# plot_grid(vp1, vp2, vp3, ncol = 3)
