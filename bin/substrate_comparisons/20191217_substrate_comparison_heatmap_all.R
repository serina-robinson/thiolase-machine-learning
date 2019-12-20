# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", "randomcoloR", 
               "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "calculated_slopes"
fils

ll <- fils[grepl(suffix, fils)]
ll <- ll[!grepl("BocPhe|furf|scratch|only|averaged|benzoate|round|C6|2019-11-15_7Ph", ll)] # rep1|rep2|reps|round|
ll

rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(cmpnd =  paste0(word(word(filename, 2, sep = "output\\/\\/"), 1, sep = "\\/"), "_",
                                word(word(filename, 2, sep = "_"), 1, sep = "_")))
  #dplyr::mutate(cmpnd = case_when(grepl("C6", filename) ~ "hexanoate",
   #                               TRUE ~ word(word(filename, 2, sep = "_"), 1, sep = "_")))

# Write to file
maprdat_long <- rawdat %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope))
summary(maprdat_long$log_slope)

# Find the average
maprdat_avg <- maprdat_long %>%
  dplyr::group_by(cmpnd) %>% # this would necessarily be (org, cmpnd) once we have biological replicates!!!!
  dplyr::summarise_each(funs(mean), log_slope)

summary(maprdat_avg$log_slope)

ggplot(aes(x = cmpnd, y = log_slope), data = maprdat_long) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(width=.1, height=0)) +
  theme_bw()

# write_csv(maprdat_avg, "output/all_cmpnds_avg_log_slopes.csv")

maprdat_merg <- as.data.frame(maprdat_long, stringsAsFactors = F)

# Convert to wide format
maprdat_wide <- reshape2::dcast(maprdat_merg, org ~ cmpnd, value.var = "log_slope") 
maprdat_wide[is.na(maprdat_wide)] <- 0

rawdat_mat <- maprdat_wide %>%
  dplyr::select(-org) %>%
  # dplyr::select(dat_order) %>%
  as.matrix()

# Set color palette
pal <- inferno(80)
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
rownames(dedup_sort)[grep("Mycolici", rownames(dedup_sort))]
newnames <- lapply(
  rownames(dedup_sort),
  function(x) bquote(italic(.(x))))

pdf("output/substrate_comparisons/substrate_comparison_heatmap_unclustered_log10_scale_per_hr_no_cutoff_sorted.pdf", width = 11, height = 11)
pheatmap(
  cluster_cols = T,
  cluster_rows = T,
  border_color = NA,
  mat = dedup_sort, 
  color = pal2,
  annotation_names_row = T, 
  fontsize_col = 10, 
  fontsize_row = 8, 
  annotation_names_col =  T,
  labels_row = as.expression(newnames))
dev.off()

#colnames(dedup)[grep("C6", colnames(dedup))] <- c("hexanoate averaged")
pdf("output/substrate_comparisons/substrate_comparison_heatmap_unclustered_transposed_log10_scale_per_hr_no_cutoff.pdf", width = 20, height = 8)
pheatmap(
  cluster_cols = F,
 # cluster_rows = F,
  border_color = NA,
  mat = t(dedup_sort), 
  #  breaks = mat_breaks,
  color = pal2,
  annotation_names_row = T, 
  fontsize_col = 8, 
  fontsize_row = 10, 
  annotation_names_col = T,
  labels_col = as.expression(newnames))
dev.off()



