# Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer", 
               "ggplot2", "data.table", "maditr", "broom")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

### 2019-11-08
date <- "2019-11-08"
cmpnd <- "heptynoate"

# Read in the heptynoate slopes biorep 1
hept1 <- read_csv("output/2019-10-25/2019-10-25_heptynoate_all_data_calculated_slopes.csv")

# Read in the heptynoate slopes biorep 2 
hept2 <- read_csv('output/2019-11-08/2019-11-08_heptynoate_biorep2_all_data_calculated_slopes.csv')

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "calculated_slopes"
fils

# Select the files of interest
ll0 <- fils[grepl(suffix, fils)]
ll <- ll0[grepl("heptynoate", ll0)]
ll <- ll[!grepl("_only_", ll)]

# Read in the data
rawdat <- tibble(filename = ll) %>%
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(cmpnd = case_when(TRUE ~ paste0(word(word(filename, 2, sep = "_"), 1, sep = "_"), "_", word(filename, 3, sep = "_"))))
rawdat$cmpnd <- gsub("_all", "_biorep1", rawdat$cmpnd)

# Calculate log slope
maprdat_long <- rawdat %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope)) %>%
  dplyr::select(cmpnd, org, log_slope)
# summary(maprdat_long$log_slope)
maprdat_long

# Now combine with all organisms not present/active at all
orgkey <- read_csv("data/72_OleA_masterwell_org_key.csv")
orgkey$orgs <- gsub("_", " ", orgkey$orgs)
orgkey$orgs
inactives <- orgkey$orgs[!orgkey$orgs %in% maprdat_long$org[maprdat_long$cmpnd == "heptynoate_biorep1"]]
inactives
inactive_df1 <- data.frame(cmpnd = rep("heptynoate_biorep1", length(inactives)), inactives, rep(0, length(inactives)), stringsAsFactors = F)
colnames(inactive_df1) <- c("cmpnd", "org", "log_slope")

inactives <- orgkey$orgs[!orgkey$orgs %in% maprdat_long$org[maprdat_long$cmpnd == "heptynoate_biorep2"]]
inactive_df2 <- data.frame(cmpnd = rep("heptynoate_biorep2", length(inactives)), inactives, rep(0, length(inactives)), stringsAsFactors = F)
colnames(inactive_df2) <- c("cmpnd", "org", "log_slope")

 #colnames = colnames(maprdat_mat))

# Find the average
maprdat_avg <- maprdat_long %>%
  bind_rows(inactive_df1) %>%
  bind_rows(inactive_df2) %>%
  dplyr::group_by(org) %>%
  dplyr::summarise_each(funs(mean), log_slope) %>%
  dplyr::mutate(cmpnd = "heptynoate_average") %>%
  dplyr::select(cmpnd, org, log_slope)

# To write to file
maprdat_write <- maprdat_avg %>%
  dplyr::filter(log_slope != 0) %>%
  write_csv(., "output/scratch_output/heptynoate_averaged_calculated_slopes.csv")

maprdat_join <- maprdat_long %>%
  bind_rows(maprdat_avg)
 #  dplyr::left_join(., maprdat_avg, by = "org")

maprdat_merg <- as.data.frame(maprdat_join, stringsAsFactors = F)

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

dedup_sort <- dedup[order(rowSums(dedup), decreasing = T),]
head(dedup_sort)

pdf("output/substrate_comparisons/substrate_comparison_heatmap_unclustered_heptynoate_only_log10_scale_per_hr_no_cutoff_sorted.pdf", width = 3.5, height = 9)
pheatmap(
  cluster_cols = F,
  cluster_rows = F,
  border_color = NA,
  mat = dedup_sort, 
  color = pal2,
  annotation_names_row = T, 
  fontsize_col = 10, 
  fontsize_row = 8, 
  annotation_names_col =  T)
dev.off()


pdf("output/substrate_comparisons/substrate_comparison_heatmap_unclustered_heptynoate_only_transposed_log10_scale_per_hr_no_cutoff.pdf", width = 10, height = 3)
pheatmap(
  cluster_cols = F,
  cluster_rows = F,
  border_color = NA,
  mat = t(dedup_sort), 
  #  breaks = mat_breaks,
  color = pal2,
  annotation_names_row = T, 
  fontsize_col = 8, 
  fontsize_row = 10, 
  annotation_names_col =  T)
dev.off()

