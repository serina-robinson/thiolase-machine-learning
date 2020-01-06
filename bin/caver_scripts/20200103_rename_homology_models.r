# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest","gbm","nnet","e1071","svmpath","lars",
               "glmnet","svmpath", "data.table", "bio3d")


# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the summary info about the homology models
summ <- fread("data/homology_models/summaryinfo_extended", sep = "|", fill = T, data.table = F, header = T) %>%
  janitor::clean_names() %>%
  dplyr::filter(!grepl("Resolution", resolution)) %>%
  dplyr::group_by(number_job_description) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(org = paste0(word(number_job_description, sep = "_1_", 2))) %>%
  dplyr::mutate(pdb_id = gsub("PDBTitle: ", "", v12)) %>%
  dplyr::mutate(pdb_id = gsub("crystal structure of ", "", pdb_id)) %>%
  dplyr::mutate(filnams = paste0(job_id, ".final.pdb"))
  #dplyr::filter(pdb_id == "3-oxoacyl-(acyl carrier protein)2 synthase iii, fabh (xoo4209) from xanthomonas oryzae pv.3 oryzae kacc10331")



# Read in the old homology models
fils <- data.frame(list.files("data/homology_models/", pattern = ".pdb"), stringsAsFactors = F)
colnames(fils) <- "filnams"

comb$number_job_description

comb <- summ %>%
  dplyr::left_join(fils, ., by = "filnams") %>%
  dplyr::mutate(new_filnam = paste0(number_job_description, ".pdb")) %>%
  dplyr::filter(pdb_id == "3-oxoacyl-(acyl carrier protein)2 synthase iii, fabh (xoo4209) from xanthomonas oryzae pv.3 oryzae kacc10331")


for(i in 1:nrow(comb)) {
  pdb_temp <- read.pdb(paste0("data/homology_models/", comb$filnams[i]))
  write.pdb(pdb_temp, file = paste0("data/caver_models/homology_models/", comb$new_filnam[i]))
}

# Read in the new homology models
seventytwo <- list.files("data/caver_models/homology_models/")
write_csv(data.frame(seventytwo, stringsAsFactors = F), "data/caver_models/CastP_submission_log.csv")
# newls <- list.files("data/13_homology_models_X_oryzae_template/")

