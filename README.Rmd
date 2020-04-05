---
title: "Example analysis"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library("tidyverse")
library("viridis")
library("readxl")
library("ggtree")
library("caret")
library("rsample")
library("ranger")
library("pROC")
library("RColorBrewer")
library("ggpubr")
library("ggpmisc")
library("Biostrings")
library("DECIPHER")
library("kableExtra")
```

## 

Welcome to the Github repo for scripts to reproduce figures and analysis in the following manuscript: 

Robinson, S.L., Smith, M.D., Richman, J.E., Aukema, K.G., & Wackett, L.P. (2020) Machine learning-based prediction of activity and substrate specificity for OleA enzymes in the thiolase superfamily. **Under review.**

Below are some example scripts of analysis conducted in this study

### 1. Reading in and cleaning up the data 

```{r}
maprdat <- read_csv("data/enzyme_substrate_activity_data.csv")
head(maprdat)

# Find the average across biological triplicates
maprdat_avg <- maprdat %>%
  dplyr::group_by(cmpnd, org) %>% 
  dplyr::summarise_each(list(mean_activity = mean), log_slope) 
```

### 2. Visualize phylogenetic tree of enzyme sequences paired with enzyme activity heatmap

```{r, echo = F}
# Convert to wide format
maprdat_merg <- as.data.frame(maprdat_avg, stringsAsFactors = F)
maprdat_wide <- reshape2::dcast(maprdat_merg, org ~ cmpnd, value.var = "mean_activity") 
maprdat_wide[is.na(maprdat_wide)] <- 0

rawdat_mat <- maprdat_wide %>%
  dplyr::select(-org) %>%
  as.matrix()

# Set color palette
pal <- inferno(80)
pal2 <- pal[c(10:80)]

# Fix names
maprdat_mat <- rawdat_mat
rownames(maprdat_mat) <- maprdat_wide$org
rownames(maprdat_mat) <- gsub("_", " ", rownames(maprdat_mat))
rownames(maprdat_mat) <- gsub("XC", "Xanthomonas campestris OleA*", rownames(maprdat_mat))
rownames(maprdat_mat) <- gsub("Pseudoxanthomonas", "Pseudoxanthomonas sp.", rownames(maprdat_mat))
rownames(maprdat_mat) <- paste0(word(rownames(maprdat_mat), 1, sep = " "), " ", word(rownames(maprdat_mat), 2, sep = " "))

# Remove duplicates and exceptions
dedup <- maprdat_mat[!duplicated(rownames(maprdat_mat)),]
dedup_sort <- dedup[order(rowSums(dedup), decreasing = T),]

# Add in benzoate
testr <- cbind(dedup_sort, rep(0, nrow(dedup_sort)))
colnames(testr)[ncol(testr)] <- "benzoate"

testm <- testr
testm[testm > 0] <- 1
nsub <- data.frame(cbind(rownames(testm), rowSums(testm)))
colnames(nsub) <- c("org", "nsub")

tr <- "data/73_OleA_aligned_trimmed_fasttree_phylo.nwk"
fastree <- treeio::read.tree(tr)
ggtree_phylo <- treeio::as.phylo(fastree)

# Read in the taxonomy
rawtax <- read_excel("data/OleA_taxonomic_classification.xlsx") 
xantho <- rawtax[grep("Xanthomonas", rawtax$organism),]
xantho$organism <- c("Xanthomonas_campestris_OleA")
xantho$genus <- "Xanthomonas_campestris"

tax <- rawtax  %>%
  bind_rows(xantho)
tax$organism <- gsub(" ", "_", tax$organism)
ptax <- ggtree_phylo$tip.label

for(i in 1:length(tax$organism)) {
  ind <- grep(tax$organism[i], ptax)
  tax$acc[i] <- word(ptax[ind], sep = "\\.1", 1)
  ptax[ind] <- paste0(ptax[ind], "_", tax$class[i], "_", tax$phylum[i])
}
tax$acc[is.na(tax$acc)] <- "NP_635607"

# Reorder the data frame to match
dt2 <- data.frame(ptax, word(ptax, sep = "\\.1", 1))
colnames(dt2) <- c("label", "acc")

# Merge with the tax df
merg <- dt2 %>%
  inner_join(tax, by = "acc") %>%
  dplyr::mutate(labs = paste0(word(organism, sep = "_", 1), " ",
                word(organism, sep = "_", 2)))

# Fix discrepancies in merging names
pseudo1 <- rownames(maprdat_mat)[grep("Pseudoxanthomonas", rownames(maprdat_mat))]
pseudo2 <- merg$labs[grep("Pseudoxanthomonas", merg$labs)][1]
merg$labs[grep("Pseudoxanthomonas", merg$labs)] <- pseudo1

leif1 <- rownames(maprdat_mat)[grep("Leifsonia", rownames(maprdat_mat))]
leif2 <- merg$labs[grep("Leifsonia", merg$labs)][1]
merg$labs[grep("Leifsonia", merg$labs)] <- leif1

# Create a phylogenetic tree
ggtree_phylo$tip.label <- merg$labs
namord <- names(sort(colMeans(testr), decreasing = T))
test_ord <- testr[,match(namord, colnames(testr))]

ggt <- ggtree(ggtree_phylo) + 
  geom_tiplab(aes(fontface = "italic"),
              align = TRUE, linetype = 2, offset = 0.04) + 
  xlim(NA, 20)

# Calculate average activity
avg_act <- test_ord %>%
  rowMeans(.) %>%
  data.frame(., stringsAsFactors = F) %>%
  rownames_to_column(., var = "org")
colnames(avg_act) <- c("org", "avg_activity")

# Merge with taxonomic information
mergtax <- avg_act %>%
  left_join(., merg, by = c("org" = "labs"))

# Read in the custom palette
clustkey <- read_csv("data/OleA_palette_key.csv")
clustkey$levs[clustkey$levs == "Green non-sulfur bacteria"] <- "Chloroflexi"

# Fix discrepancies in naming
mergtax$class[grepl("Opit", mergtax$class)] <- "Opitutae"
mergord <- mergtax[rev(match(ggt$data$label[ggt$data$isTip], mergtax$org)),] %>%
  dplyr::select(-label) %>%
  dplyr::mutate(label = org) %>%
  dplyr::left_join(., clustkey, by = c("class" = "levs")) %>%
  dplyr::select(label, org, avg_activity, acc, organism, genus, family, order, class, phylum, pal2)
df2 <- data.frame(cbind(mergord$acc, mergord$organism, mergord$label, as.numeric(mergord$avg_activity)), stringsAsFactors = F)

df3 <- df2[order(df2$X2, decreasing = T),] %>%
  dplyr::left_join(., nsub, by = c("X3" = "org")) %>%
  dplyr::mutate(organism = gsub("_", " ", X2)) %>%
  dplyr::mutate(activity = round(as.numeric(X4), 3)) %>%
  dplyr::select(X1, organism, activity, nsub) %>%
  arrange(desc(activity))
colnames(df3) <- c('NCBI Accession', 'Organism', 'Average activity', 'Total number of substrates accepted')

ggdf <- ggt %<+% mergord
clustord <- clustkey[match(levels(as.factor(ggdf$data$class)), clustkey$levs),]

# Draw a phylogenetic tree
gsz <- ggdf +
  geom_tippoint(aes(size = avg_activity), color = ggdf$data$pal2[ggdf$data$isTip], 
                fill = ggdf$data$pal2[ggdf$data$isTip], x = 10.35)
```



```{r a_taller_figure, fig.height = 20, fig.width = 11}
# Pair phylogenetic tree with heatmap
gheatmap(gsz, test_ord, offset = 5.75, width = 2, 
         colnames_position = "top", colnames_offset_y = 2,
         colnames_angle = 60, font.size = 4, color = NULL) +
  scale_fill_viridis(option = "inferno") +
  theme(legend.title = element_blank(), 
        legend.position = "none")
```

```{r legend, fig.height = 3, fig.width = 2.5, echo = F}
# Taxonomic legend
plot.new()
par(mar=c(0.01,0.01,0.01,0.01))
legend("center", bty = "n", bg = "white",
       legend = clustord$levs,  pch = 19, 
       cex = 1, pt.cex = 1, 
       col = clustord$pal2)
```

### 3. Machine learning classification of active vs. inactive enzyme-substrate pairs
```{r}
# Read in the features for machine learning 
dat <- read_csv("data/machine_learning/20200111_1095_training_examples_12angstrom_prot_features.csv")

# Split into test and training sets
dat_split <- initial_split(dat, strata = "is_active")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 % training, 25 % testing

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("id", "is_active")]
x_test <- dat_test[,!colnames(dat_test) %in% c("id", "is_active")]
y_train <- dat_train$is_active
y_test <- dat_test$is_active

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)

# Random forest 10-fold cross-validation (not run due to computational time)
# rf <- train(
#   x = df_train,
#   y = y_train,
#   method = "ranger",
#   metric = "ROC",
#   trControl = trainControl(method = "repeatedcv", number = 10,
#                            repeats = 3,
#                            verboseIter = T, classProbs = T,
#                            savePredictions = "final",
#                            returnResamp = "all"),
#   num.trees = 1000,
#   verbose = TRUE,
#   importance = "permutation") 
# saveRDS(rf, "data/machine_learning/models/20200111_rf_10foldcv.rds")

# Read in the random forest model from cross-validation
rf_bin_xval <- readRDS("data/machine_learning/models/20200111_rf_10foldcv.rds")
rf_bin_xval$bestTune # Hyperparameters from cross-validation

# Receiver operating characteristic curve
rfRoc <- roc(response = rf_bin_xval$pred$obs,
             predictor = rf_bin_xval$pred$Y,
             levels = rev(levels(rf_bin_xval$pred$obs)))
par(bty = "L")
plot(rfRoc, type = "s",
     col = "#529DCF", xaxs = "i", yaxs="i",
     print.auc = TRUE, print.auc.x = 0.8, print.auc.y = 0.6)
```
```{r, echo = F}
vdf <- varImp(rf_bin_xval, scale = T)
dtfvec <- vdf$importance[order(vdf$importance, decreasing = T),]
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

reds <- rev(c("firebrick4", "firebrick1", "lightpink"))
blues <- rev(c("midnightblue", "royalblue1", "lightskyblue1"))
blues <- brewer.pal(9, "Blues")[c(3, 6, 9)]

vpavg <- ggplot(data = vimp_combined,
                aes(x=reorder(var_fix, importance), y=importance, fill=color_fct))+
  geom_bar(stat="identity", position ="dodge",  color = "black") + 
  ylab("Variable Importance")+
  xlab("")+
  guides(fill = F) +
  scale_fill_manual(values = blues) +
  theme_pubr()+
  coord_flip()

vpavg
```

### 4. Machine learning regression of enzyme-substrate activity levels

```{r}
# Set random seed 
set.seed(1234)

# Read in the activity data
dat <- read_csv('data/machine_learning/20200111_550_regression_training_examples_12angstrom_prot_features.csv')

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

# Random forest regression 10-fold cross-validation (not run due to computational time)
# rf_mod <- train(
#   x = df_train,
#   y = y_train,
#   method = "ranger",
#   trControl = trainControl(method = "repeatedcv", number = 10,
#                            repeats = 3,
#                            savePredictions = "final"),
#   num.trees = 1000,
#   verbose = TRUE,
#   importance = "permutation")
# saveRDS(rf_mod, "data/machine_learning/models/220200104_final_rf_regression_model.rds")

# rf_reg_xval <- readRDS("data/machine_learning/models/20200104_final_rf_regression_model.rds")
rf_reg_xval <- readRDS("data/machine_learning/models/20200111_rf_mod_10foldcv.rds")

# Calculate root-mean square error
Metrics::rmse(rf_reg_xval$pred$pred, rf_reg_xval$pred$obs) 

# rf_reg <- ranger(y_train ~., data = form_train, num.trees = 1000, 
#                  splitrule = rf_reg_xval$bestTune$splitrule,
#                  mtry = rf_reg_xval$bestTune$mtry, 
#                  min.node.size = rf_reg_xval$bestTune$min.node.size,
#                  importance = "permutation")
# 
# rf_pred <- predict(rf_reg, data = form_test)$predictions
# Metrics::rmse(y_test, rf_pred) # 0.219
# rf_df <- data.frame(cbind(rf_pred, y_test))

my.formula <- y ~ x
summary(lm(rf_reg_xval$pred$pred ~ rf_reg_xval$pred$obs))

rf_df <- data.frame(cbind(rf_reg_xval$pred$pred, rf_reg_xval$pred$obs))
rf_df_resid <- rf_df %>%
  mutate(resid = rf_reg_xval$pred$pred - rf_reg_xval$pred$obs)

colnames(rf_df_resid) <- c("pred", "obs", "resid")

ggplot(rf_df_resid, aes(x = pred, y = resid)) +
  geom_point(alpha = .3) + 
  geom_line(col = "red", y = 0.0,
              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Predicted enzyme activity") +
  ylab("Residuals") 

ggplot(rf_df_resid, aes(x = obs, y = pred)) +
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", method = "lm",
              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Observed enzyme activity") +
  ylab("Predicted enzyme activity") +
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..rr.label.., sep = "~~~")),
               parse = TRUE)
```

```{r, echo = F}
rf_reg_xval <- readRDS("data/machine_learning/models/20200111_rf_mod_10foldcv.rds")

vdf <- varImp(rf_reg_xval, scale = T)
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
  theme_pubr() +
  coord_flip()

vpavg
```

### 5. Prediction for a new sequence

```{r}
# Align active site residues and extract
new_query_seq <- "data/new_test_OleA.fasta"
nqs <- readAAStringSet("data/new_test_OleA.fasta")
nqs$Vulcaniibacterium_new_seq 

source("lib/extract_12angstrom_aas.R")
extract_84_list <- lapply(1:length(length(new_query_seq)), 
                          function(x) { extract_12angstrom_aas(new_query_seq) })
extract_84_df <- data.frame(matrix(unlist(extract_84_list), nrow = length(extract_84_list), byrow=T), 
                stringsAsFactors = FALSE)

# Set column names
colnams84 <- read_csv("shiny_app/data/seqft_colnames.csv", col_names = F) %>%
  pull()
colnames(extract_84_df) <- colnams84

# Remove variables with non-zero variance
which_rem <- read_csv("shiny_app/data/vars_to_remove.csv", col_names = T) %>%
  pull()

# Combine with substrate features
molec_fts <- read_csv("shiny_app/data/PC7_molecular_descriptors_clean.csv")
hex <- molec_fts[,2:8] 
subs <- rep(pull(molec_fts[,1]), length(new_query_seq))
hex_df <- data.frame(cbind(subs, apply(hex, 2, function(x) rep(x, length(new_query_seq)))), stringsAsFactors = F)
colnames(hex_df) <- c("substrate", paste0("PC", 1:7))


hex_df_sorted <- hex_df[order(hex_df$PC1),] 
ex84_df_scaled <- data.frame(apply(extract_84_df, 2, function(x) rep(x, 15)), stringsAsFactors = F)

test_df <- bind_cols(hex_df_sorted, ex84_df_scaled) %>%
dplyr::select(-all_of(which_rem), -substrate)
test_df[,1:7] <- apply(test_df[,1:7], 2, as.numeric)

# Binary classification
rf_funct_class <- readRDS("shiny_app/data/models/20200324_rf_bin_class_all_data_no_prot_unscaled_optimized.rds")
rf_funct_class_pred <- predict(rf_funct_class, data = test_df, predict.all = F)
rf_fc_pred <- rf_funct_class_pred$predictions
res_fc_prob <- apply(rf_fc_pred, 1, max)
res_fc <- tibble(colnames(rf_fc_pred)[apply(rf_fc_pred, 1, which.max)])

# Return a data frame of predictions
pred_df <- data.frame(rep(names(nqs), 15), hex_df_sorted$substrate, res_fc, round(res_fc_prob, 2), stringsAsFactors = F) 

# Regression
rf_reg <- readRDS("shiny_app/data/models/20200401_rf_regression_all_data_no_prot_unscaled_optimized.rds")
rf_reg <- predict(rf_reg, data = test_df, predict.all = F)
rf_reg_pred <- rf_reg$predictions

reg_df <- data.frame(cbind(pred_df, round(rf_reg_pred, 2)), stringsAsFactors = F)
colnames(reg_df) <- c("query_name", "subst", "is_active", "pred_prob", "raw_activity")
```

```{r, echo = F}
# Show picture of substrate
reg_df$substrate <- paste0('<img src = "https://github.com/serina-robinson/thiolase-machine-learning/raw/master/shiny_app/www/', ... = reg_df$subst, '.jpeg" style="width: 90%;  height: auto;" >')

reg_df$raw_activity[reg_df$is_active == "N"] <- "No predicted activity"  
ord_df <- reg_df[order(reg_df$query_name),] %>%
  dplyr::mutate(IUPAC_name = case_when(subst == "decanoate" ~ "4-nitrophenyl decanoate",
                    subst == "dodecanoate" ~ "4-nitrophenyl dodecanoate",
                    subst == "7Ph.heptanoate" ~ "4-nitrophenyl 7-phenylheptanoate",
                    subst == "hexanoate" ~ "4-nitrophenyl hexanoate",
                    subst == "heptanoate" ~ "4-nitrophenyl heptanoate",
                    subst == "biotin" ~ "4-nitrophenyl biotin",
                    subst == "heptynoate" ~ "4-nitrophenyl heptynoate",
                    subst == "oxadiazole" ~ "4-nitrophenyl 3-(5-phenyl-1,3,4-oxadiazol-2-yl)propanoate",
                    subst == "dimethyl" ~ "4-nitrophenyl 2,2-dimethylhexanoate",
                    subst == "butoxy" ~ "4-nitrophenyl 2-(2-butoxyethoxy)acetate",
                    subst == "cyclopentyl" ~ "4-nitrophenyl 3-cyclopentylpropanoate",
                    subst == "TMA" ~ "4-nitrophenyl trimethylacetate",
                    subst == "azido" ~ "4-nitrophenyl 6-azidohexanoate",
                    subst == "benzoate" ~ "4-nitrophenyl benzoate",
                    subst == "ClPh.propionate" ~ "4-nitrophenyl-3-(4-chlorophenoxy)propanoate",
                    TRUE ~ "")) %>%
dplyr::select(query_name, IUPAC_name, substrate, is_active,	pred_prob, raw_activity)

ord_df %>%
    mutate(`Query name` = query_name) %>%
    mutate(`Substrate` = substrate) %>%    
    mutate(`IUPAC name` = IUPAC_name) %>%
    mutate(`Active?` = cell_spec(is_active, color = "white",
                               background = case_when(is_active == "Y" & pred_prob >= 0.65 ~ "green",
                                                      is_active == "N" & pred_prob >= 0.65 ~ "darkred",
                                                      pred_prob < 0.65 ~ "grey"))) %>%
    mutate(`Prediction probability` = cell_spec(pred_prob, color = "white",
                               background = case_when(is_active == "Y" & pred_prob >= 0.65 ~ "green",
                                                      is_active == "N" & pred_prob >= 0.65 ~ "darkred",
                                                      pred_prob < 0.65 ~ "grey"))) %>%
    mutate(`Estimated log10 activity ☨` = cell_spec(raw_activity, color = "white",
                                background = case_when(is_active == "Y" & pred_prob >= 0.65 ~ "green",
                                                    is_active == "N" & pred_prob >= 0.65 ~ "darkred",
                                                    pred_prob < 0.65 ~ "grey"))) %>%
    select(-query_name, -pred_prob, -is_active, -raw_activity, -substrate, -IUPAC_name) %>%
    kable(., format = 'html', escape = F) %>%
    kable_styling(font_size = 15,
        bootstrap_options = c("striped", "hover", "condensed"))
```