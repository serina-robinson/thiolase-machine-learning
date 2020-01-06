# Read in the raw data
# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "gbm","nnet","e1071","svmpath","lars",
               "glmnet", "svmpath", "readxl", "ggpubr", "plotROC", "ROCR")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the principal componenets of molecular features
molec_fts <- read_csv("data/machine_learning/PC7_molecular_descriptors.csv") %>%
  dplyr::rename(substrate = V1)
molec_fts$substrate <- gsub(" ", "\\.", molec_fts$substrate)
molec_fts$substrate[molec_fts$substrate == "oxidazole"] <- "oxadiazole"

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
  dplyr::select(id, contains("PC"), contains("4KU5"), is_active)

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
                 splitrule = "extratrees",
                 mtry = 279, min.node.size = 1,
                 importance = "permutation")
  
  rf_models[[i]]$models <- rf_1
  rf_models[[i]]$oob <- 1 - rf_1$prediction.error # Classifcation accuracy 1 - OOB error
  
  rf_1_pred <- predict(rf_1, form_test)

  rf_models[[i]]$confmat <- confusionMatrix(rf_1_pred$predictions, as.factor(dat_test$is_active))
  
  vimp_1 <- data.frame(cbind(sort(rf_1$variable.importance, decreasing = T),
                             names(sort(rf_1$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
    dplyr::rename(importance = X1,
                  variable = X2) %>%
    dplyr::mutate(importance = as.numeric(importance)) %>%
    dplyr::slice(1:30)

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
mean_testing_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$confmat$overall[1]})
mean_training_vec
mean_testing_vec

conftab <- data.frame(list(rf_models[[which.max(mean_testing_vec)]]$confmat$table))

pdf("output/machine_learning/best_rf_model_confmat.pdf", width = 3, height = 3)
ggplot(conftab, aes(Prediction, Reference)) +
  geom_tile(aes(fill = Freq)) + 
  geom_text(aes(label = round(Freq, 1))) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank()) +
  theme_pubr()
dev.off()
rf_models[[which.max(mean_testing_vec)]] # 82.42 % testing accuracy

saveRDS(rf_models[[which.max(mean_testing_vec)]]$models, "data/machine_learning/models/final_random_forest_classification_model.rds")


dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) 
sd_df <- dat_df %>%
  dplyr::summarise_each(sd)
sd_df
colMeans(dat_df)


rf_test_pred <- predict(rf_models[[which.max(mean_testing_vec)]]$models, form_test)

comp_df <- data.frame(cbind(dat_test$id, rf_test_pred$predictions, y_test)) %>%
  dplyr::mutate(y_pred = as.numeric(V2)) %>%
  dplyr::mutate(y_test = as.numeric(y_test)) %>%
  dplyr::mutate(incorrect = abs(y_pred - y_test)) %>%
  dplyr::select(V1, y_test, y_pred, incorrect) %>%
  dplyr::filter(incorrect == 1)
comp_df



### VIMP
vimp_combined <- lapply(1:length(rf_models), function(x) { rf_models[[x]]$vimp }) %>%
  map_dfr(~ .) %>%
  group_by(variable) %>%
  summarise_at(vars(importance), mean) %>%
  arrange(desc(importance)) %>%
  dplyr::mutate(var_fix = case_when(grepl("PC3", variable) ~ "Aromaticity (PC3)",
                                    grepl("PC1", variable) ~ "Molecular weight (PC1)",
                                    grepl("PC2", variable) ~ "Molecular topology (PC2)",
                                    grepl("PC4", variable) ~ "Solubility (PC4)",
                                    grepl("PC5", variable) ~ "Oxygen content (PC5)",
                                    grepl("PC6", variable) ~ "Nitrogen content (PC6)",
                                    grepl("PC7", variable) ~ "Chlorine content (PC7)",
                TRUE ~ paste0(word(variable, -1, sep = "_"), " ", word(variable, 1, sep = "_"))))
vimp_combined



vpavg <- ggplot(data = vimp_combined[1:30,],
              aes(x=reorder(var_fix, importance), y=importance, fill=importance))+
  geom_bar(stat="identity", position="dodge") + coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F) +
  theme_pubr()
vpavg  

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
                 splitrule = "extratrees",
                 mtry = 279, min.node.size = 1,
                 importance = "permutation",
                 probability = T)
  
  rf_models[[i]]$models <- rf_1
  rf_models[[i]]$oob <- 1 - rf_1$prediction.error # Classifcation accuracy 1 - OOB error
  
  rf_1_pred <- predict(rf_1, form_test)
  
  finpreds <- as.factor(ifelse(rf_1_pred$predictions[,2] > 0.5, "Y", "N"))
  
  rf_models[[i]]$confmat <- confusionMatrix(finpreds, as.factor(dat_test$is_active))
  
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
mean_testing_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$confmat$overall[1]})
mean_training_vec
mean_testing_vec

conftab <- data.frame(list(rf_models[[which.max(mean_testing_vec)]]$confmat$table))


ggplot(conftab, aes(Prediction, Reference)) +
  geom_tile(aes(fill = Freq)) + 
  geom_text(aes(label = round(Freq, 1))) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank()) +
  theme_pubr()


dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) 
sd_df <- dat_df %>%
  dplyr::summarise_each(sd)
sd_df
colMeans(dat_df)

rf_mod_best <- rf_models[[which.max(mean_testing_vec)]]$models
df_all <- data.frame(cbind(dat_train$id, rf_mod_best$predictions, y_train)) %>%
  dplyr::mutate(y_pred = as.numeric(V2)) %>%
  dplyr::mutate(y_train = as.numeric(y_train)) %>%
  dplyr::mutate(incorrect = y_pred - y_train) %>%
  dplyr::mutate(false_pos = incorrect == 1) %>%
  dplyr::mutate(false_neg = incorrect == -1)
writeLines(as.character(df_all$V1[df_all$false_neg == T]))
writeLines(as.character(df_all$V1[df_all$false_pos == T]))

wridf_all$V1[df_all$false_neg == T]


rf_test_pred <- predict(rf_models[[which.max(mean_testing_vec)]]$models, form_test)
finpreds <- as.factor(ifelse(rf_test_pred$predictions[,2] > 0.5, "Y", "N"))
data.frame(cbind(dat_test$id, finpreds, y_test))

comp_df <- data.frame(cbind(dat_test$id, rf_test_pred$predictions[,2], finpreds, y_test)) %>%
  dplyr::mutate(y_pred = as.numeric(finpreds)) %>%
  dplyr::mutate(y_test = as.numeric(y_test)) %>%
  #dplyr::rename(y_pred = V2) %>%
  dplyr::mutate(incorrect = abs(y_pred - y_test)) %>%
  dplyr::mutate(y_test_bool = case_when(y_test == 2 ~ 1,
                              y_test == 1 ~ 0)) %>%
  dplyr::rename(y_pred_prob = V2) %>%
  dplyr::select(V1, y_pred_prob, y_test_bool) 
  #dplyr::select(V1, y_test, y_pred, incorrect) %>%
  # dplyr::filter(incorrect == 1) 



writeLines(as.character(comp_df$V1))

data(ROCR.simple)


comp_df$y_pred_prob <- as.numeric(comp_df$y_pred_prob)


pdf("output/machine_learning/random_forest_ROC_curve.pdf", width = 5, height = 5)
rocplot <- ggplot(comp_df, aes(m = y_pred_prob, d = y_test_bool)) + 
  geom_roc(labels=FALSE, n.cuts=0) +
  theme_pubr() 
rocplot +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC =", round(calc_auc(rocplot)$AUC, 3))) +
  xlab("False positive fraction") +
  ylab("True positive fraction")
dev.off()
calc_auc(rocplot)

# rocplot + style_roc(theme = theme_grey) #+ #geom_rocci(fill="pink") 



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
