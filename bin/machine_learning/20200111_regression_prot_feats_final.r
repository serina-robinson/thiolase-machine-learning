# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "gbm","nnet","e1071","svmpath","lars", "ggpmisc",
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

# Set random seed 
set.seed(1234)

# Split into test and training data
dat_split <- initial_split(dat, strata = "activity")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("id", "activity")]
x_test <- dat_test[,!colnames(dat_test) %in% c("id", "activity")]
y_train <- dat_train$activity
y_test <- dat_test$activity

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)

y_train

# Random forest regression
rf_mod <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           savePredictions = "final"),
  # tuneGrid = tgrid,
  num.trees = 1000,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  importance = "permutation") # BEST RMSE is 55.7%

# Confusion matrix
getTrainPerf(rf_mod) # RMSE is 25%
rf_imp <- varImp(rf_mod)

ggplot(rf_imp, top = 32) + xlab("")
rf_mod$bestTune

# Predict for new data
rf_pred <- predict(rf_mod, newdata = form_test)
Metrics::rmse(y_test, rf_pred) # 0.219
rf_df <- data.frame(cbind(rf_pred, y_test))

my.formula <- y ~ x
summary(lm(rf_df$y_test ~ rf_pred))
summary(lm(rf_pred ~ rf_df$y_test))


rf_df_resid <- rf_df %>%
  mutate(resid = y_test - rf_pred)

pdf("output/machine_learning/20200111_random_forest_regression_testing_results_residual_plot.pdf", width = 4, height = 4)
ggplot(rf_df_resid, aes(x = rf_pred, y = resid)) +
  # geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  geom_line(col = "red", y = 0.0,
              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Predicted enzyme activity") +
  ylab("Residuals") 
  # stat_poly_eq(formula = my.formula, 
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #              parse = TRUE) 

dev.off()

rf_pred
pdf("output/machine_learning/20200111_random_forest_regression_testing_results_observed_predicted.pdf", width = 4, height = 4)
ggplot(rf_df, aes(x = y_test, y = rf_pred)) +
  # geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", method = "lm",
              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Observed enzyme activity") +
  ylab("Predicted enzyme activity") +
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..rr.label.., sep = "~~~")),
               parse = TRUE)

dev.off()

saveRDS(rf_mod, "data/machine_learning/models/20200111_rf_regression_10foldcv.rds")

rfdf_train <- data.frame(cbind(rf_mod$finalModel$predictions, y_train), stringsAsFactors = F)
colnames(rfdf_train) <- c("Pred", "Obs")

pdf("output/machine_learning/20200111_random_forest_regression_training_results_observed_predicted.pdf", width = 4, height = 4)
ggplot(rfdf_train, aes(x = Pred, y = Obs)) +
  # geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", method = "lm",
              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Observed enzyme activity") +
  ylab("Predicted enzyme activity") +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE) 
dev.off()

rf_mod$finalModel$variable.importance
rf_mod$finalModel$variable.importance




# dtf <- data.frame(cbind(names(rf_mod$finalModel$variable.importance), rf_mod$finalModel$variable.importance), stringsAsFactors = F)
# dtfsort <- dtf[order(dtf$X2, decreasing = T),]
# dtfsort[1:10,]
vdf <- varImp(rf_mod, scale = T)
dtfvec <-vdf$importance[order(vdf$importance, decreasing = T),]
dtfnams <- rownames(vdf$importance)[order(vdf$importance, decreasing = T)]
dtfsort <- data.frame(cbind(dtfnams[1:30], dtfvec[1:30]), stringsAsFactors = F)

colnames(dtfsort)[1:2] <- c("variable", "importance")
dtfsort$variable <- as.character(dtfsort$variable)
dtfsort$importance <- as.numeric(dtfsort$importance)

vimp_combined <- dtfsort %>%
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
                                    grepl("PP3", variable) ~ "Protein H-bonding (Cruciani property 3)",
                                    grepl("Z1", variable) ~ "Protein lipophilicity index",
                                    grepl("KF6", variable) ~ "Protein partial specific volume",
                                    grepl("ST3", variable) ~ "Structural-topology scale factor 3",
                                    grepl("F5", variable) ~ "Local flexibility (FASGAI factor 5)",
                                    grepl("pi_prot", variable) ~ "Protein isoelectric point")) %>%
  dplyr::mutate(color_fct = case_when(grepl("PC", var_fix) ~ "3",
                                      grepl("\\(aa", var_fix) ~ "2",
                                      TRUE ~ "1")) 

  
# grepl("^KF", variable) ~ gsub("KF", "Kidera factor ")))
vimp_combined$variable[grepl("secstr", vimp_combined$variable)]
vimp_combined

reds <- rev(c("firebrick4", "firebrick1", "lightpink"))
blues <- rev(c("midnightblue", "royalblue1", "lightskyblue1"))
blues <- brewer.pal(9, "Blues")[c(3, 6, 9)]


vpavg <- ggplot(data = vimp_combined,
                aes(x=reorder(var_fix, importance), y=importance, fill=color_fct))+
  geom_bar(stat="identity", position="dodge",  color = "black") + 
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F) +
  scale_fill_manual(values = reds) +
  theme_pubr()+
#  scale_y_continuous( limits = c(0, 0.08)) +
  coord_flip()
#coord_cartesian(xlim=c(0,0.08))
vpavg

pdf("output/machine_learning/20200113_rf_regression_caret_train_varimp.pdf", height = 7, width = 5)
vpavg
dev.off()

vimp30 <- vimp_combined[1:30,]
vimp30
imps <- vimp30$var_fix[grepl("aa", vimp30$var_fix)]
imps <- word(imps, sep = "aa ", -1)
imps <- gsub("\\)", "", imps)
unique(imps)
paste0(unique(imps), collapse = "+")

# D151, # ALA 173

# Now test elastic net (glmnet)
enet_mod <- train(
  x = df_train,
  y = y_train,
  method = "glmnet",
  trControl = trainControl(method = "repeatedcv", 
                           repeats = 3, 
                           number = 10,
                           savePredictions = "final"))

enet_mod$results # best RMSE is 0.276 which is pretty good!

# Look at important coefficients
important <- names(enet_mod$finalModel$coefficients)[!is.na(enet_mod$finalModel$coefficients)]
which_coeffs <- important[2:length(important)]
enet_imp <- varImp(enet_mod)
ggplot(enet_imp, top = 50) + xlab("")
getTrainPerf(enet_mod)

# TrainRMSE TrainRsquared  TrainMAE method
# 0.2746507     0.5491163 0.2136555 glmnet
saveRDS(enet_mod, "data/machine_learning/models/20200111_enet_mod_10foldcv.rds")

# Predict with MARS model for new data
enet_pred <- predict(enet_mod, newdata = form_test)
Metrics::rmse(y_test, enet_pred) #0.278
enet_df <- data.frame(cbind(enet_pred, y_test))

my.formula <- y ~ x
summary(lm(enet_df$y_test ~ enet_pred))
pdf("output/machine_learning/20200111_enet_regression_testing_results_observed_predicted.pdf", width = 5, height = 5)
ggplot(rf_df, aes(x = y_test, y = enet_pred)) +
  # geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", method = "lm",
              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Observed enzyme activity") +
  ylab("Predicted enzyme activity") +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) 

dev.off()


## Multivariate adaptive regression splines
mars_grid <- expand.grid(degree = 1:2, nprune = seq(2, 26, by = 2))

mars_mod <- train(
  x = df_train,
  y = y_train,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = trainControl(method = "repeatedcv", 
                           repeats = 3, 
                           number = 10,
                           savePredictions = "final"))

# Training and test performance
getTrainPerf(mars_mod) # RMSE is 0.284
mars_mod$finalModel
mars_mod$bestTune
summary(mars_mod$finalModel)
saveRDS(mars_mod, "data/machine_learning/models/20200111_mars_mod_10foldcv.rds")
mars_pred <- predict(mars_mod, newdata = form_test)
Metrics::rmse(y_test, mars_pred)
mars_mod$finalModel$coefficients

# Variable importance
mars_imp <- varImp(mars_mod)
ggplot(mars_imp, top = 10) + xlab("")
channel_a <- c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353)
channel_b <- c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246)
channels <- c(channel_a, channel_b)
imps <- word(rownames(mars_imp$importance)[mars_imp$importance != 0], -1, sep = "_")
imps <- as.numeric(imps[!grepl("PC", imps)])
imps[imps %in% channel_a]
imps[imps %in% channel_b]


# Predict for new data
mars_pred <- predict(mars_mod, newdata = form_test) 
Metrics::rmse(y_test, mars_pred) #0.252
mars_df <- data.frame(cbind(mars_pred, y_test))

my.formula <- y ~ x
summary(lm(mars_df$y_test ~ mars_pred))
pdf("output/machine_learning/20200111_mars_regression_testing_results_observed_predicted.pdf", width = 5, height = 5)
ggplot(rf_df, aes(x = y_test, y = mars_pred)) +
  # geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", method = "lm",
              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Observed enzyme activity") +
  ylab("Predicted enzyme activity") +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) 

dev.off()

