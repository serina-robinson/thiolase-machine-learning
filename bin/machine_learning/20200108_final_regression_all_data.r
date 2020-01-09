# Read in the raw data
# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "gbm","nnet","e1071","svmpath","lars",
               "glmnet", "svmpath", "readxl", "ggpubr", "ggpmisc")

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

p1 <- ggplot(dedup) +
  geom_histogram(aes( x = activity, y = ..density..), alpha = 0.7,
                 binwidth = 0.075, fill = "gray80", color = "black") +
  theme_pubr(base_size = 17) +
  xlab("Enzyme activity") +
  ylab("Density")

pdf("output/machine_learning/full_activity_distribution.pdf")
p1
dev.off()

# write_csv(dedup, "data/machine_learning/20191228_1095_training_examples_12angstrom_features.csv")
dat <- dedup %>%
  dplyr::mutate(id = paste0(org, "_", substrate)) %>%
  dplyr::mutate(is_active = as.factor(case_when(activity == 0 ~ "N",
                                                TRUE ~ "Y"))) %>%
  dplyr::filter(is_active == "Y") %>%
  dplyr::select(-which_rem, -org, -substrate) %>%
  dplyr::select(id, contains("PC"), contains("4KU5"), activity)



# Set random seed 
set.seed(1234)

x_all <- dat %>%
  dplyr::select(-id, -activity)
y_all <- dat$activity
form_all <- data.frame(cbind(x_all, y_all), stringsAsFactors = F, row.names = dat$id)
head(form_all)

rf_1 <- ranger(y_all ~., data = form_all, num.trees = 1000, 
               splitrule = "extratrees",
               mtry = 136, min.node.size = 5,
               importance = "permutation")

rmse <- sqrt(rf_1$prediction.error) # Classifcation accuracy 1 - OOB error

df2 <- data.frame(cbind(rf_1$predictions, y_all))
colnames(df2) <- c("Predicted", "Truth")

my.formula = y ~ x
pdf("output/machine_learning/obs_pred_regression_plot_for_fig3_all.pdf", width = 4, height = 4)
ggplot(df2, aes(x = Truth, y = Predicted)) +
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", method = "lm", 
              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr(base_size = 14) #+
  # stat_poly_eq(formula = my.formula,
  #              aes(label = paste(..rr.label..)),
  #              parse = TRUE)
dev.off() 
