## Example machine learning workflow to predict enzyme substrate specificity

Welcome to an example workflow* to generate figures and analysis in the following manuscript: 

Robinson, S.L., Smith, M.D., Richman, J.E., Aukema, K.G., & Wackett, L.P. (2020) Machine learning-based prediction of activity and substrate specificity for OleA enzymes in the thiolase superfamily. **Under review.**

&ast; **Note:** this is an abbreviated version of the workflow. For a complete RMarkdown version of this document go to [https://rpubs.com/robinsonserina/thiolases](https://rpubs.com/robinsonserina/thiolases)

```{r}
# Install and load packages
pacman::p_load("tidyverse", viridis", "readxl", "ggtree", "caret", "rsample",
"ranger", "pROC", "RColorBrewer", "ggpubr", "ggpmisc", "Biostrings", "DECIPHER","kableExtra")
```

### 1. Reading in and cleaning up the data 

```{r}
maprdat <- read_csv("data/enzyme_substrate_activity_data.csv")

# Find the average across biological triplicates
maprdat_avg <- maprdat %>%
  dplyr::group_by(cmpnd, org) %>% 
  dplyr::summarise_each(list(mean_activity = mean), log_slope) 
```

### 2. Visualize phylogenetic tree of enzyme sequences paired with enzyme activity heatmap

```{r}
# Pair phylogenetic tree with heatmap
gheatmap(gsz, test_ord, offset = 5.75, width = 2, 
         colnames_position = "top", colnames_offset_y = 2,
         colnames_angle = 60, font.size = 4, color = NULL) +
  scale_fill_viridis(option = "inferno") +
  theme(legend.title = element_blank(), 
        legend.position = "none")
```
<img src="https://github.com/serina-robinson/thiolase-machine-learning/raw/master/pngs/fig2.png">

##### Taxonomic legend
<img src="https://github.com/serina-robinson/thiolase-machine-learning/raw/master/pngs/tax_legend.png" height="200" width="250">

### 3. Machine learning classification of active vs. inactive enzyme-substrate pairs

```{r}
# Read in the features for machine learning 
dat <- read_csv("data/machine_learning/20200111_1095_training_examples_12angstrom_prot_features.csv")

# Set seed 
set.seed(1234)

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

# Random forest 10-fold cross-validation
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
<img src="https://github.com/serina-robinson/thiolase-machine-learning/raw/master/pngs/auroc.png" height="400" width="600">

##### Variable importance

<img src="https://github.com/serina-robinson/thiolase-machine-learning/raw/master/pngs/vimp1.png" height="500" width="600">

### 4. Machine learning regression of enzyme-substrate activity levels

```{r}
# Set seed 
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

# Random forest regression 10-fold cross-validation 
rf_mod <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           savePredictions = "final"),
  num.trees = 1000,
  verbose = TRUE,
  importance = "permutation")
# saveRDS(rf_mod, "data/machine_learning/models/220200104_final_rf_regression_model.rds")

# rf_reg_xval <- readRDS("data/machine_learning/models/20200104_final_rf_regression_model.rds")
rf_reg_xval <- readRDS("data/machine_learning/models/20200111_rf_mod_10foldcv.rds")

# Calculate training set root-mean square error 
Metrics::rmse(rf_reg_xval$pred$pred, rf_reg_xval$pred$obs) 

rf_reg <- ranger(y_train ~., data = form_train, num.trees = 1000,
                 splitrule = rf_reg_xval$bestTune$splitrule,
                 mtry = rf_reg_xval$bestTune$mtry,
                 min.node.size = rf_reg_xval$bestTune$min.node.size,
                 importance = "permutation")

rf_pred <- predict(rf_reg, data = form_test)$predictions
Metrics::rmse(y_test, rf_pred) # 0.219
rf_df <- data.frame(cbind(rf_pred, y_test))

my.formula <- y ~ x
rf_df_resid <- rf_df %>%
  mutate(resid = y_test - rf_pred)
colnames(rf_df_resid) <- c("pred", "obs", "resid")

ggplot(rf_df_resid, aes(x = pred, y = resid)) +
  geom_point(alpha = .3) +
  geom_line(col = "red", y = 0.0,
              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Predicted enzyme activity") +
  ylab("Residuals")
```
<img src="https://github.com/serina-robinson/thiolase-machine-learning/raw/master/pngs/resid.png" height="400" width="600">

```
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
<img src="https://github.com/serina-robinson/thiolase-machine-learning/raw/master/pngs/resid2.png" height="400" width="600">

##### Variable importance
<img src="https://github.com/serina-robinson/thiolase-machine-learning/raw/master/pngs/vimp2.png"
height="500" width="600">

### 5. Prediction of substrate specificity for a new protein sequence
```{r}

# Align with template sequence and extract active site residues
new_query_seq <- "data/new_test_OleA.fasta"
nqs <- readAAStringSet("data/new_test_OleA.fasta")

source("lib/extract_12angstrom_aas.R")
extract_84_list <- lapply(1:length(length(new_query_seq)), 
                          function(x) { extract_12angstrom_aas(new_query_seq) })
extract_84_df <- data.frame(matrix(unlist(extract_84_list), nrow = length(extract_84_list), byrow=T), 
                stringsAsFactors = FALSE)
colnams84 <- read_csv("shiny_app/data/seqft_colnames.csv", col_names = F) %>%
  pull()
colnames(extract_84_df) <- colnams84

# Remove sequence features with non-zero variance
which_rem <- read_csv("shiny_app/data/vars_to_remove.csv", col_names = T) %>%
  pull()

# Create a data frame of substrate features
molec_fts <- read_csv("shiny_app/data/PC7_molecular_descriptors_clean.csv")
hex <- molec_fts[,2:8] 
subs <- rep(pull(molec_fts[,1]), length(new_query_seq))
hex_df <- data.frame(cbind(subs, apply(hex, 2, function(x) rep(x, length(new_query_seq)))), stringsAsFactors = F)
colnames(hex_df) <- c("substrate", paste0("PC", 1:7))
hex_df_sorted <- hex_df[order(hex_df$PC1),] 
ex84_df_scaled <- data.frame(apply(extract_84_df, 2, function(x) rep(x, 15)), stringsAsFactors = F)

# Create full data frame of features
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

<img src="https://github.com/serina-robinson/thiolase-machine-learning/raw/master/pngs/tbl.png" alt="heatmap">