# Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", 
               "RColorBrewer", "svglite", "ggplot2",
               "data.table", "pheatmap", "gplots", "DECIPHER")

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
  dplyr::mutate(pdb_id = gsub("crystal structure of ", "", pdb_id))
unique(summ$pdb_id)

summ$pdb_id

non_x <- summ %>%
  dplyr::filter(pdb_id != "3-oxoacyl-(acyl carrier protein)2 synthase iii, fabh (xoo4209) from xanthomonas oryzae pv.3 oryzae kacc10331")
non_x$org

# Find the sequences
sqs <- readAAStringSet("data/homology_models/72_OleA_JGI_unaligned_for_Phyre.fasta")

orgs <- paste0(word(names(sqs), sep = "_", 4), "_", word(names(sqs), sep = "_", 5))
orgs[grep("Bacillus_sp", orgs)] <- "Bacillus_sp__7586_K"
# orgs[grep("Bacillus", orgs)] <- Bacillus_sp__7586_K

which_temps <- sqs[orgs %in% non_x$org]
writeXStringSet(which_temps, "data/homology_models/13_OleA_with_diff_templates_for_Phyre.fasta")
                
                
mymat <- as.matrix(table(summ$pdb_id, summ$org)) #59 have the following 


rownames(mymat)

pdf("output/72_homology_model_template_heatmap.pdf", width = 15, height = 40)
#par(mar = c(1, 1, 1, 1))
heatmap(t(mymat), keep.dendro = FALSE, margins = c(30, 15), cexRow = 0.75, labCol = c("X. oryzae OleA",
                                                                            "Burkholderia xenovorans unknown",
                                                                            "Streptomyces peucetius PKS",
                                                                            "Propionibacterium Kas III", 
                                                                            "MxnB (M.fulvus)", 
                                                                            "S. aureus Kas III"
                                                                            ))
dev.off()
hmp
# ggsave(hmp, "output/heatmap.pdf", device = "svg")
