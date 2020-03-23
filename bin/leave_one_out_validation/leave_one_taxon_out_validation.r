# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", "readxl", "ggpubr", "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the data
dat <- read_csv("data/machine_learning/20200111_1095_training_examples_12angstrom_prot_features.csv") %>%
  dplyr::mutate(spec = word(id, sep = "_", 1))
dat$spec

# Read in the taxonomic identification
rawtax <- read_excel("data/OleA_taxonomic_classification.xlsx") 
xantho <- rawtax[grep("Xanthomonas", rawtax$organism),]
xantho$organism <- c("Xanthomonas_campestris_OleA")
xantho$genus <- "Xanthomonas_campestris"

tax <- rawtax  %>%
  bind_rows(xantho)

tax$organism <- gsub("_", " ", tax$organism)
dtax <- tax %>%
  dplyr::mutate(spec = paste0(word(organism, sep = " ", 1), " ", word(organism, sep = " ", 2)))

# Fix discrepancies in merging names
pseudo1 <- dtax$spec[grep("Pseudoxanthomonas", dtax$spec)]
pseudo2 <- dat$spec[grep("Pseudoxanthomonas", dat$spec)] 
dat$spec[grep("Pseudoxanthomonas", dat$spec)] <- pseudo1

leif1 <- dtax$spec[grep("Leifsonia", dtax$spec)]
leif2 <- dat$spec[grep("Leifsonia", dat$spec)] 
dat$spec[grep("Leifsonia", dat$spec)] <- leif1

# Merge with the tax df
merg <- dtax %>%
  inner_join(dat, by = "spec") %>%
  dplyr::select(-organism, -genus, -family, -order, -phylum, -spec)
uniq_tax <- unique(merg$class)

# Set random seed 
set.seed(5678)  

# Create a list of models
rf_models <- vector(mode = "list",
                    length = length(uniq_tax))
uniq_tax

# Leave-one-compound-out validation
for(i in 1:length(uniq_tax)) {
  
  # Identify compound to be 'left out' of model training
  loo_tax <- uniq_tax[i]
  
  # Leave one compound out
  dat_train <- merg[merg$class != loo_tax,]
  dat_test  <- merg[merg$class == loo_tax,]
  print(dim(dat_test))
  
  # Define variables of interest
  x_train <- dat_train %>%
    dplyr::select(-id, -is_active, -class)
  x_test <- dat_test %>%
    dplyr::select(-id, -is_active, -class)
  
  y_train <- dat_train$is_active
  y_test <- dat_test$is_active
  
  df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)
  
  ## Train model using default parameters (gini split-rule, sqrt(m) for mtry)
  rf_1 <- ranger(y_train ~., data = form_train, num.trees = 1000, 
                 splitrule = "gini",
                 mtry = sqrt(ncol(form_train)),
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

# Benzoate causes an error because the model has not been trained on that type of chemical features previously

# Will report the results for the 14 other compounds
mean_training_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$oob }) 
mean_training_vec

mean_testing_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$confmat$overall[1]})
accuracy <- round(mean_testing_vec * 100, 2)
accuracy

dat_df <- tibble(uniq_tax) %>%
  bind_cols(tibble(accuracy)) %>%
  arrange(desc(accuracy))
dat_df

colnames(dat_df) <- c("left-out taxon", "leave-one-taxon-out accuracy (%)")
write_csv(dat_df, "output/Leave-one-taxon-out-analysis.csv")

writeLines(uniq_tax)




