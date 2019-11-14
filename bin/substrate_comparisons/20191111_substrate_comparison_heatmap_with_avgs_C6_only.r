# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", "randomcoloR", 
               "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "calculated_slopes"

# Select the files of interst
ll0 <- fils[grepl(suffix, fils)]
ll <- ll0[grepl("hexanoate|C6", ll0)]
ll <- ll[!grepl("avged", ll)]
ll
rawdat <- tibble(filename = ll) %>%
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(date_run = substr(word(filename, 1, sep = "_"), nchar("output//2019-11-08/") + 1, nchar("output//2019-11-08/") + 10)) %>%
  dplyr::mutate(cmpnd = case_when(grepl("round1", filename) ~ "hexanoate round 1",
                                  grepl("round2", filename) ~ "hexanoate round 2",
                                  TRUE ~ "hexanoate round 3"))

# Write to file
maprdat_log <- rawdat %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope)) 
summary(maprdat_long$log_slope)

# Find the average
maprdat_avg <- maprdat_log %>%
  dplyr::group_by(org) %>%
  dplyr::summarise_each(funs(mean), log_slope) #
maprdat_avg$cmpnd <- "hexanote averaged"

# write_csv(maprdat_avg, "output/C6_avg_log_slopes.csv")
# Combine with original
maprdat_long <- maprdat_log %>%
 bind_rows(., maprdat_avg)
 
maprdat_merg <- as.data.frame(maprdat_long, stringsAsFactors = F)

# Convert to wide format
maprdat_wide <- reshape2::dcast(maprdat_merg, org ~ cmpnd, value.var = "log_slope") 
maprdat_wide[is.na(maprdat_wide)] <- 0

rawdat_mat <- maprdat_wide %>%
  dplyr::select(-org) %>%
  # dplyr::select(dat_order) %>%
  as.matrix()

# Set color palette
pal <- magma(80)
pal2 <- pal[c(10:80)]

# Fix names
maprdat_mat <- rawdat_mat
rownames(maprdat_mat) <- maprdat_wide$org
rownames(maprdat_mat) <- gsub("_", " ", rownames(maprdat_mat))
rownames(maprdat_mat) <- gsub("XC", "Xanthomonas campestris OleA*", rownames(maprdat_mat))
rownames(maprdat_mat) <- gsub("Pseudoxanthomonas", "Pseudoxanthomonas sp.", rownames(maprdat_mat))
rownames(maprdat_mat) <- paste0(word(rownames(maprdat_mat), 1, sep = " "), " ", word(rownames(maprdat_mat), 2, sep = " "))
head(maprdat_mat)

# Now combine with all organisms not present/active at all
orgkey <- read_csv("data/72_OleA_masterwell_org_key.csv")
inactives <- orgkey$orgs[!orgkey$orgs %in% rownames(maprdat_mat)]
inactive_mat <- matrix(ncol = ncol(maprdat_mat), nrow = length(inactives), 0) #colnames = colnames(maprdat_mat))
colnames(inactive_mat) <- colnames(maprdat_mat)
rownames(inactive_mat) <- inactives

full_mat <- rbind(maprdat_mat, inactive_mat)
rownames(full_mat) <- gsub("_", " ", rownames(full_mat))
rownames(full_mat) <- paste0(word(rownames(full_mat), 1, sep = " "), " ", word(rownames(full_mat), 2, sep = " "))

# Remove duplicates and exceptions
dedup <- full_mat[!duplicated(rownames(full_mat)),]
# dedup_df <- data.frame(dedup, stringsAsFactors = F)
dedup <- dedup[rownames(dedup) != "Pseudoxanthomonas NA",]
dedup <- dedup[rownames(dedup) != "Lysobacter tolerans",]

# dedup_sort <- dedup[order(dedup[,4], decreasing = T),]
head(dedup_sort)
dedup_sort <- dedup[order(rowSums(dedup), decreasing = T),]

pdf("output/substrate_comparisons/substrate_comparison_heatmap_unclustered_C6_only_log10_scale_per_hr_no_cutoff_sorted.pdf", width = 3.5, height = 9)
pheatmap(
  cluster_cols = F,
  cluster_rows = F,
  border_color = NA,
  mat = dedup_sort[,1:3], 
  color = pal2,
  annotation_names_row = T, 
  fontsize_col = 10, 
  fontsize_row = 8, 
  annotation_names_col =  T)
dev.off()


pdf("output/substrate_comparisons/substrate_comparison_heatmap_unclustered_C6_only_transposed_log10_scale_per_hr_no_cutoff.pdf", width = 10, height = 3)
pheatmap(
  cluster_cols = F,
  cluster_rows = F,
  border_color = NA,
  mat = t(dedup_sort[,1:3]), 
  #  breaks = mat_breaks,
  color = pal2,
  annotation_names_row = T, 
  fontsize_col = 8, 
  fontsize_row = 10, 
  annotation_names_col =  T)
dev.off()



