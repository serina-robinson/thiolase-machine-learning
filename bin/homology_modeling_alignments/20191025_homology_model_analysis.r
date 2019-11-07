# Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer", "svglite", "ggplot2", "data.table", "pheatmap", "gplots")

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
summ$org

mymat <- as.matrix(table(summ$pdb_id, summ$org)) #59 have the following 
unique(summ$pdb_id)

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
