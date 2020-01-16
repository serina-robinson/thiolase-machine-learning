# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "PerformanceAnalytics",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", "corrr", "plyr",
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
  dplyr::select(-which_rem, -enzyme, -is_active, -sqs, -core, -acc, -nams) #-IUPAC, -SMILES, -cmpnd_abbrev, ) 

# Center and scale everything except for id and is_active
datscale <- dat %>%
  select(-id, -activity, -org, -substrate) %>%
  scale(., center = TRUE) %>%
  data.frame(., stringsAsFactors = F) %>%
  mutate_if(is.character,as.numeric)

dat <- datscale %>%
  bind_cols(tibble(dat$id), tibble(dat$org), tibble(dat$substrate), .) %>%
  bind_cols(., tibble(dat$activity))

colnames(dat)[1:3] <- c("id", "org", "substrate")
colnames(dat)[ncol(dat)] <- "activity"

whichmax <- dat %>%
  dplyr::group_by(substrate) %>%
  dplyr::arrange(desc(activity)) %>%
  dplyr::slice(1:3) %>%
  dplyr::select(org, substrate, activity) %>%
  dplyr::arrange(org) %>%
  dplyr::filter(grepl("Kytococcus", org))
whichmax

pdf("output/machine_learning/KF_correlations.pdf", height = 20, width = 2)
ggplot(dat, aes(x = KF1, y = activity)) + 
  geom_point() + 
  geom_smooth(method = "loess", colour = "red", fill = "red") + 
  geom_smooth(method = "lm", colour = "blue", fill = "blue") + 
  facet_grid(as.factor(dat$substrate)) +
  # (data=cors, aes(label=paste("r=", cor, sep="")), x=1, y=1) +
  theme_pubr()


cor.test(dat$KF1, dat$activity)
cor.test(dat$hydrophob_core, dat$activity)

cors <- ddply(dat, c("substrate"), summarise, cor = round(cor(KF1, activity), 2))
cors[order(cors$cor, decreasing = T),]

dat$hydrophob_core
dat$activity
cor_hyd <- ddply(dat, c("substrate"), summarise, cor = round(hydrophob_core, activity), 2)

pdf("output/machine_learning/hydrophob_core_correlations.pdf", height = 20, width = 2)
ggplot(dat, aes(x = hydrophob_core, y = activity)) + 
  geom_point() + 
  geom_smooth(method = "loess", colour = "red", fill = "red") + 
  geom_smooth(method = "lm", colour = "blue", fill = "blue") + 
  facet_grid(as.factor(dat$substrate)) +
  # (data=cors, aes(label=paste("r=", cor, sep="")), x=1, y=1) +
  theme_pubr()
dev.off()

p + facet_grid(vs ~ am) +
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=30, y=4)

avg_dat <- dat %>%
  # dplyr::select(org, substrate, KF1, activity) %>%
  group_by(org) %>%
  # dplyr::select(-org, -substrate)
  dplyr::summarise_each(funs(mean), activity) %>%
  dplyr::rename(avg_activity = activity)

avg_dat

avg_sub_org_dat <- dat %>%
  group_by(org, substrate) %>%
  dplyr::summarise_each(funs(mean), activity) %>%
  dplyr::rename(avg_sub_activity = activity)

avgdf2 <- avg_sub_org_dat %>%
  dplyr::left_join(., dat, by = "org") %>%
  dplyr::select()

head(avg_dat)

datavg_comb <- avg_dat %>%
  dplyr::left_join(., dat, by = "org") %>%
  #dplyr::select(-contains("PC"), -contains("4KU5"), -id, -substrate, -org, -activity) %>%
  distinct() %>%
  dplyr::select(contains("PC"), activity)
  #::select(KF1, F2, ST1, hmoment, mw_prot, lngth_prot, boman_interactions, hydrophob_core, avg_activity)
datavg_comb
chart.Correlation(datavg_comb)


cor.test(avg_comb$avg_activity, avg_comb$KF1)



# MOre alpha helix and bend properties
pdf("output/substrate_comparisons/activity_with_structural_correleations.pdf")
chart.Correlation(avg_comb, pch = 19, histogram = TRUE)
dev.off()
