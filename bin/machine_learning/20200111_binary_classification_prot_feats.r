# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "gbm","nnet","e1071","svmpath","lars", "pROC",
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
colnames(seq_fts)

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
  dplyr::select(-which_rem, -org, -substrate, -enzyme, -activity, -sqs, -core, -acc, -nams) #-IUPAC, -SMILES, -cmpnd_abbrev, ) 

# Center and scale everything except for id and is_active
datscale <- dat %>%
  select(-id, -is_active) %>%
  scale(., center = TRUE) %>%
  data.frame(., stringsAsFactors = F) %>%
  mutate_if(is.character,as.numeric)

dat <- datscale %>%
  bind_cols(tibble(dat$id), .) %>%
  bind_cols(., tibble(dat$is_active))
colnames(dat)[1] <- "id"
colnames(dat)[ncol(dat)] <- "is_active"

# Set random seed 
set.seed(1234)

# Split into test and training data
table(dat$is_active)
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

# Random forest
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  metric = "ROC",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final",
                           returnResamp = "all"),
  num.trees = 1000,
  verbose = TRUE,
  importance = "permutation") 

# Confusion matrix
getTrainPerf(rf) # Training set accuracy is 82% for random forest
# saveRDS(rf, "data/machine_learning/models/20200111_rf_10foldcv.rds")
rf$pred$pred
rf$pred$obs
rf$pred$Y

rfRoc <- roc(response = rf$pred$obs,
             predictor = rf$pred$Y,
             levels = rev(levels(rf$pred$obs)))
rfRoc

pdf("output/machine_learning/rocPlot.pdf",
    width=4, height=4)
par(bty = "L")
plot(rfRoc, type = "s", 
     col = "#529DCF", xaxs = "i", yaxs="i",
     print.auc = TRUE, print.auc.x = 0.8, print.auc.y = 0.6)
dev.off()

sink("output/machine_learning/20200111_rf_binary_classification_results.txt")
cm_rf
sink()

vimp
vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T),
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:15)
vimp

#pdf("output/machine_learning/random_forest_12angstroms_activity_binary_classification_vimp.pdf", width = 6, height = 6)
ggplot(data = vimp,
       aes(x=reorder(variable,importance), y=importance, fill=importance))+
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
#dev.off()

### Random forest results
getTrainPerf(rf)


mod_train_pred <- as.factor(ifelse(rf$finalModel$predictions[,1] >= 0.5, "N", "Y"))
train_cm <- confusionMatrix(mod_train_pred, y_train)
cmdf <- data.frame(train_cm$table) %>%
  dplyr::mutate(Truth = ifelse(Reference == "N", "inactive", "active")) %>%
  dplyr::mutate(`Model Prediction` = ifelse(Prediction == "N", "inactive", "active"))
cmdf

pdf("output/machine_learning/rf_training_set_confMat.pdf", width = 3, height = 2.5)
ggplot(cmdf, aes(`Model Prediction`, Truth)) +
  geom_tile(aes(fill = Freq)) + 
  geom_text(aes(label = round(Freq, 1))) +
  scale_fill_gradient(low = "white", high = "#4393C7") +
  theme_pubr() +
  theme(legend.position = "none") 
dev.off()

rf$bestTune$splitrule
rf$bestTune$mtry
rf$bestTune$min.node.size

rf_ml <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = as.character(rf$bestTune$splitrule),
                mtry = rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
                importance = "permutation", probability = T)
rf$bestTune
rf_pred <- predict(rf, newdata = form_test)

cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$is_active))
cm_rf # 83.88% test set classification accuracy

cm_rf_test <- data.frame(cm_rf$table) %>%
  dplyr::rename(Truth = Reference) %>%
  dplyr::rename(`Model Prediction` = Prediction)

pdf("output/machine_learning/rf_testing_set_confMat.pdf", width = 2.5, height = 2.5)
ggplot(cm_rf_test, aes(`Model Prediction`, Truth)) +
  geom_tile(aes(fill = Freq)) + 
  geom_text(aes(label = round(Freq, 1))) +
  scale_fill_gradient(low = "white", high = "#4393C7") +
  theme_pubr() +
  theme(legend.position = "none") 
dev.off()

vimp <- data.frame(cbind(sort(rf$finalModel$variable.importance, decreasing = T),
                         names(sort(rf$finalModel$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:40)
ggplot(data = vimp,
       aes(x=reorder(variable,importance), y=importance, fill=importance))+
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue") +
  theme_pubr()


vdf <- varImp(rf, scale = T)
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
                                    grepl("T4", variable) ~ "Topological descriptor 4",
                                    grepl("F3", variable) ~ "Bulky properties (FASGAI factor 3)",
                                    grepl("ST6", variable) ~ "Structural-topology scale factor 6")) %>%
  dplyr::mutate(color_fct = case_when(grepl("PC", var_fix) ~ "3",
                                      grepl("\\(aa", var_fix) ~ "2",
                                      TRUE ~ "1")) 


# grepl("^KF", variable) ~ gsub("KF", "Kidera factor ")))
vimp_combined$variable[grepl("secstr", vimp_combined$variable)]
vimp_combined

reds <- rev(c("firebrick4", "firebrick1", "lightpink"))
blues <- rev(c("midnightblue", "royalblue1", "lightskyblue1"))
blues <- brewer.pal(9, "Blues")[c(3, 6, 9)]

pdf("output/machine_learning/20200113_rf_classification_caret_train_varimp.pdf", height = 7, width = 5)
vpavg <- ggplot(data = vimp_combined,
                aes(x=reorder(var_fix, importance), y=importance, fill=color_fct))+
  geom_bar(stat="identity", position="dodge",  color = "black") + 
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F) +
  scale_fill_manual(values = blues) +
  theme_pubr()+
  #  scale_y_continuous( limits = c(0, 0.08)) +
  coord_flip()
#coord_cartesian(xlim=c(0,0.08))
vpavg
dev.off()

vimp30 <- vimp_combined[1:30,]
vimp30
imps <- vimp30$var_fix[grepl("aa", vimp30$var_fix)]
imps <- word(imps, sep = "aa ", -1)
imps <- gsub("\\)", "", imps)
unique(imps)
paste0(unique(imps), collapse = "+")


pdf("output/machine_learning/20200113_rf_regression_caret_train_varimp.pdf", height = 7, width = 5)
vpavg
dev.off()

#pdf("output/machine_learning/random_forest_12angstroms_activity_binary_classification_vimp.pdf", width = 6, height = 6)



###Naive Bayes method
#nb_grid <- expand.grid(usekernel = TRUE, fL = 0, adjust = 1)

nb <- train(
  x = df_train, 
  y = y_train,
  method = "nb",
  #tuneGrid = nb_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, 
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"))

# Confusion matrix
getTrainPerf(nb) # Only 59% training accuracy with NB
nb$bestTune
nb$results

nb_pred <- predict(nb, newdata = form_test)
nb$bestTune
cm_nb <- confusionMatrix(nb_pred, as.factor(dat_test$is_active))
cm_nb

sink("output/machine_learning/20200111_nb_activity_binary_classification_results.txt")
cm_nb
sink()

saveRDS(nb, "data/machine_learning/models/20200111_nb_10foldcv.rds")

## Make a heatmap of confusion matrix results
cm_list <- list(cm_nb$table)
names(cm_list) <- c("nb_aa")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]

# ggsave("output/machine_learning/nb_12angstroms_activity_binary_classification_matrix_20200102.jpeg", height=7, width=7, units='in')

###Neural networks using nnet
# Currently throws an error
nnet_grid <- expand.grid(.decay = c(0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7),
                         .size = c(3, 5, 10, 20))

nnet_mod <- train(
  x = df_train,
  y = y_train,
  method = "nnet",
  MaxNWts = 1000000000,
  # size = 1,
  tuneGrid = nnet_grid,
  trControl = trainControl(method = "repeatedcv", repeats = 3,
                           number = 10,
                           verboseIter = T,
                           savePredictions = "final",
                           classProbs = TRUE))

# Confusion matrix
nnet_mod$bestTune
getTrainPerf(nnet_mod)

vimp <- data.frame(cbind(sort(nnet_mod$variable.importance, decreasing = T),
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:15)
vimp
saveRDS(nnet_mod, "data/machine_learning/models/20200111_nnet_10foldcv.rds")


# size decay
#  10 0.5

# approx_roc_curve(nnet_ml, "Neural network") %>%
#   ggplot(aes(x = 1 - specificity, y = sensitivity)) +
#   geom_path()  +
#   geom_abline(col = "red", alpha = .5)

nnet_pred <- predict(nnet_mod, newdata = form_test)
cm_nnet <- confusionMatrix(nnet_pred, as.factor(dat_test$is_active))
cm_nnet # 79% training accuracy

sink("output/machine_learning/20200111_nnet_binary_classification_results.txt")
cm_nnet
sink()
# dtl_feat_select <- data.frame(round(cm_nnet$byClass[,colnames(cm_nnet$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))

## Make a heatmap of confusion matrix results
cm_list <- list(cm_nnet$table)
names(cm_list) <- c("nnet_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("output/machine_learning/20200111_nnet_binary_classification_matrix.jpeg", height=7, width=7, units='in')



