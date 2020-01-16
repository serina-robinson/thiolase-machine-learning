# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "gbm","nnet","e1071","svmpath","lars", "ggExtra",
               "glmnet", "svmpath", "readxl", "ggpubr", "doMC", "doParallel")

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
seq_fts$enzyme[seq_fts$enzyme == "4KU5_Xanthomonas_campestris"] <- "NP_635607.1_4KU5_Xanthomonas_campestris"
seq_fts <- seq_fts %>%
  dplyr::mutate(acc = word(enzyme, sep = "\\.1", 1))

# Read in the protein features
prot_fts <- read_csv("data/machine_learning/73_overall_calculated_protein_properties.csv")
# grep(paste0(prot_fts$acc, collapse = "|"), seq_fts$enzyme)

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
  dplyr::left_join(., seq_fts, by = "org") %>%
  dplyr::left_join(., prot_fts, by = "acc")

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
  dplyr::filter(is_active == "Y") %>%
  dplyr::select(-which_rem, -org, -substrate, -enzyme, -is_active, -sqs, -core, -acc, -nams) #-IUPAC, -SMILES, -cmpnd_abbrev, ) 

# Center and scale everything except for id and is_active
datscale <- dat %>%
  select(-id, -activity) %>%
  scale(., center = TRUE) %>%
  data.frame(., stringsAsFactors = F) %>%
  mutate_if(is.character,as.numeric)

dat <- datscale %>%
  bind_cols(tibble(dat$id), .) %>%
  bind_cols(., tibble(dat$activity))
colnames(dat)[1] <- "id"
colnames(dat)[ncol(dat)] <- "activity"

# Get the best tuning parameters
# rf_20200111 <- readRDS("data/machine_learning/models/20200111_rf_regression_10foldcv.rds")
# mtry splitrule min.node.size
# 168  variance             5

# Set random seed 
# set.seed(5678) # missing A261
set.seed(1234)
# set.seed(1111)
# set.seed(2020)

rf_models <- vector(mode = "list",
                    length = 10)

# Split randomly into 10 test/training sets
for(i in 1:10) {
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
    ylab("Variable Importance") +
    xlab("") +
    guides(fill=F) +
    scale_fill_gradient(low="red", high="blue")
  rf_models[[i]]$vimp_plot <- vp1

}

mean_training_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$train_rmse }) 
mean_training_vec

mean_testing_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$test_rmse})
mean_testing_vec

rf_models[[which.min(mean_testing_vec)]] # 82.42 % testing accuracy
which.min(mean_testing_vec)
saveRDS(rf_models[[which.min(mean_testing_vec)]]$models, "data/machine_learning/models/20200104_final_rf_regression_model.rds")

dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) 
dat_df
colMeans(dat_df)  

sd_df <- dat_df %>%
  dplyr::summarise_each(sd)
sd_df

which.min(mean_testing_vec)

rf_test_pred <- predict(rf_models[[which.min(mean_testing_vec)]]$models, form_test)
rf_test_pred$predictions


df2 <- rf_test_pred$predictions %>%
  bind_cols(., y_test)


ggplot(rf_test_pred, aes(x = obs, y = pred)) +
  geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", 
              lty = 2, lwd = 1, alpha = .5)

# comp_df <- data.frame(cbind(dat_test$id, rf_test_pred$predictions, y_test)) %>%
#   dplyr::mutate(y_pred = as.numeric(V2)) %>%
#   dplyr::mutate(y_test = as.numeric(y_test))# %>%
# 
# comp_df
# #   dplyr::mutate(incorrect = abs(y_pred - y_test)) %>%
# #  # dplyr::select(V1, y_test, y_pred, incorrect) %>%
# #  # dplyr::filter(incorrect == 1)
# # comp_df
# 
# rf_models[[1]]$vimp[,1]

### VIMPS
vimp_combined <- lapply(1:length(rf_models), function(x) { rf_models[[x]]$vimp }) %>%
  map_dfr(~ .) %>%
  group_by(variable) %>%
  summarise_at(vars(importance), mean) %>%
  arrange(desc(importance)) %>%
  dplyr::mutate(var_fix = case_when(grepl("PC3", variable) ~ "Aromaticity (PC3)",
                                    grepl("PC1", variable) ~ "Molecular weight (PC1)",
                                    grepl("PC2", variable) ~ "Molecular connectivity (PC2)",
                                    grepl("PC4", variable) ~ "Solubility (PC4)",
                                    grepl("PC5", variable) ~ "Oxygen content (PC5)",
                                    grepl("PC6", variable) ~ "Nitrogen content (PC6)",
                                    grepl("PC7", variable) ~ "Chlorine content (PC7)",
                                    grepl("^KF1", variable) ~ "Helix preference (Kidera factor 1)",
                                    grepl("^F2", variable) ~ "Turn propensity (FASGAI factor 2)",
                                    grepl("hmoment", variable) ~ "Hydrophobic moment",
                                    grepl("secstr", variable) ~ paste0("Secondary structure (aa ", stringr::word(variable, -1, sep = "_"), ")"),
                                    grepl("polrty", variable) ~ paste0("Polarity (aa ", stringr::word(variable, -1, sep = "_"), ")"),
                                    grepl("molsz", variable) ~ paste0("Amino acid volume (aa ", stringr::word(variable, -1, sep = "_"), ")"),
                                    grepl("elechrg", variable) ~ paste0("Electrostatic charge (aa ", stringr::word(variable, -1, sep = "_"), ")"),
                                    grepl("^boman", variable) ~ "Protein-protein interaction index",
                                    grepl("^hmoment", variable) ~ "Hydrophobic moment",
                                    grepl("^ST1", variable) ~ "Structural-topology scale factor 1",
                                    grepl("mw_prot", variable) ~ "Protein molecular weight",
                                    grepl("lngth_prot", variable) ~ "Protein length",
                                    grepl("VHSE5", variable) ~ "Protein electronic properties",
                                    grepl("hydrophob_core", variable) ~ "Hydrophobicity (12Å protein core)",
                                    grepl("ST5", variable) ~ "Structural-topology scale factor 5",
                                    grepl("PP3", variable) ~ "Protein H-bonding (Cruciani property 3)")) %>%
  dplyr::mutate(color_fct = case_when(grepl("PC", var_fix) ~ "3",
                                      grepl("\\(aa", var_fix) ~ "2",
                                      TRUE ~ "1"))

# grepl("^KF", variable) ~ gsub("KF", "Kidera factor ")))
vimp_combined$variable[grepl("secstr", vimp_combined$variable)]
vimp_combined$variable[1:32]

reds <- rev(c("firebrick4", "firebrick1", "lightpink"))
blues <- rev(c("midnightblue", "royalblue1", "lightskyblue1"))
blues <- brewer.pal(9, "Blues")[c(3, 6, 9)]
blues


vpavg <- ggplot(data = vimp_combined[1:32,],
                aes(x=reorder(var_fix, importance), y=importance, fill=color_fct))+
  geom_bar(stat="identity", position="dodge",  color = "black") + 
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F) +
  scale_fill_manual(values = reds) +
  theme_pubr()+
  scale_y_continuous( limits = c(0, 0.08)) +
  coord_flip()
  #coord_cartesian(xlim=c(0,0.08))
vpavg
vimp_combined$importance

pdf("output/machine_learning/20200111_rf_regression_10iterations_varimp_avged.pdf", height = 7, width = 5)
vpavg
dev.off()

vimp32 <- vimp_combined[1:32,]
imps <- vimp32$var_fix[grepl("aa", vimp32$var_fix)]
imps <- word(imps, sep = "aa ", -1)
imps <- gsub("\\)", "", imps)
unique(imps)
paste0(unique(imps), collapse = ",")

channel_a <- c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353)
channel_b <- c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246)
channels <- c(channel_a, channel_b)

imps <- stringr::word(vimp_combined$variable, -1, sep = "_")
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
  dplyr::slice(1:32) %>%
  dplyr::filter(!grepl("PC", variable)) 

vimp_for_pymol
imps <- stringr::word(vimp_for_pymol$variable, -1, sep = "_")
chana <- imps[imps %in% channel_a]
chanb <- imps[imps %in% channel_b]
chans <- c(chana, chanb)
chans

no_chans <- imps[!imps %in% chans]
no_chans
paste0(unique(no_chans), collapse = "+")
paste0(unique(chans), collapse = "+")





