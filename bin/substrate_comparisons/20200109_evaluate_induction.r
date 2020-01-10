# Read in the raw data
# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "gbm","nnet","e1071","svmpath","lars",
               "glmnet", "svmpath", "readxl", "ggpubr", "plotROC", "ROCR")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")


# Read in the activity data
activity <- read_csv("data/machine_learning/20191218_all_cmpnds_avg_log_slopes_for_modeling.csv")
orgkey <- read_csv("data/72_OleA_masterwell_org_key_updated_for_Megan.csv") 
orgkey$orgs <- gsub("_", " ", orgkey$orgs)
orgkey <- orgkey %>%
  mutate(org = paste0(word(orgs, sep = " ", 1), " ", word(orgs, sep = " ", 2)))

orgkey$org <- gsub("Lysobacter", "Luteimonas", orgkey$org)
orgkey$org <- gsub("Pseudoxanthomonas NA", "Pseudoxanthomonas sp.", orgkey$org)
orgkey$org <- gsub("XC NA", "Xanthomonas campestris", orgkey$org)
orgkey$org

merg <- activity %>%
  left_join(., orgkey, by = "org")

inactives <- c("2-A7",
               "2-E10",
               "2-C8",
               "2-G10",
               "1-E3",
               "1-D2",
               "1-F5",
               "1-G2")



avg <- merg %>%
  group_by(org) %>%
  summarise_at(vars(activity), mean) %>%
  dplyr::rename(avg_activity = activity)

combavg <- merg %>%
  left_join(., avg, by = "org") %>%
  dplyr::mutate(saved = case_when(master_well %in% inactives ~ 1,
                                  TRUE ~ 0)) %>%
  write_csv(., "output/substrate_comparisons/avg_activity_all_enzymes.csv")





resdf <- merg[merg$master_well %in% inactives,]
#write_csv(resdf, "output/substrate_comparisons/avged_activities_induction_test.csv")
      