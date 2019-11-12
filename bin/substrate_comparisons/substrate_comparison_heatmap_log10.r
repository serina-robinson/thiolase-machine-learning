# Latest run: Nov 11, 2019

# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", "randomcoloR", "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "_calculated_slopes.csv"
head(fils)

ll <- fils[grepl(suffix, fils)]
ll <- ll[!grepl("BocPhe|C6|furf|rep1|rep2", ll)]
ll

rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(cmpnd = case_when(grepl("C6", filename) ~ "hexanoate",
                                  TRUE ~ word(word(filename, 2, sep = "_"), 1, sep = "_")))

# Write to file
maprdat_long <- rawdat %>%
 dplyr::filter(!grepl("output/2019-11-08", cmpnd))
maprdat_merg <- maprdat_long

# Convert to wide format
maprdat_wide <- reshape2::dcast(maprdat_merg, org ~ cmpnd, value.var = "max_slope") 
min <- min(maprdat_long$max_slope)

maprdat_wide[is.na(maprdat_wide)] <- 1e-4

rawdat_mat <- maprdat_wide %>%
  dplyr::select(-org) %>%
  # dplyr::select(dat_order) %>%
  as.matrix()

rownames(rawdat_mat) <- maprdat_wide$org

# Convert to log scale
log_10_mat <- log10(rawdat_mat)
# log_2_mat <- log(rawdat_mat, base = 2)
maprdat_mat <- log_10_mat
maprdat_mat

quantile_breaks <- function(xs, n = 20) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(maprdat_mat, n = 21)
mat_breaks
# mat_breaks <- sort(c(mat_breaks[c(1:2, 4:7)], 0.0001, 0.4, 0.6)) #0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) #0.8))#, 1.0)) #1.4))

pal <- inferno(length(mat_breaks))
pal2 <- pal[c(1, 3:length(pal))]
pal2

tax <- read_excel("data/OleA_taxonomic_classification.xlsx")
rownames(maprdat_mat)
rownames(maprdat_mat) <- gsub("_", " ", rownames(maprdat_mat))
rownames(maprdat_mat) <- gsub("XC", "Xanthomonas campestris OleA*", rownames(maprdat_mat))
rownames(maprdat_mat) <- gsub("Pseudoxanthomonas", "Pseudoxanthomonas sp.", rownames(maprdat_mat))
rownames(maprdat_mat) <- paste0(word(rownames(maprdat_mat), 1, sep = " "), " ", word(rownames(maprdat_mat), 2, sep = " "))
head(maprdat_mat)

# Now combine with all organisms not present 
orgkey <- read_csv("data/72_OleA_masterwell_org_key.csv")

inactives <- orgkey$orgs[!orgkey$orgs %in% rownames(maprdat_mat)]
inactive_mat <- matrix(ncol = ncol(maprdat_mat), nrow = length(inactives), -4) #colnames = colnames(maprdat_mat))
inactive_mat
colnames(inactive_mat) <- colnames(maprdat_mat)
rownames(inactive_mat) <- inactives

full_mat <- rbind(maprdat_mat, inactive_mat)

rownames(full_mat) <- gsub("_", " ", rownames(full_mat))
rownames(full_mat) <- paste0(word(rownames(full_mat), 1, sep = " "), " ", word(rownames(full_mat), 2, sep = " "))
dedup <- full_mat[!duplicated(rownames(full_mat)),]
dedup_df <- data.frame(dedup, stringsAsFactors = F)

dedup <- dedup[rownames(dedup) != "Pseudoxanthomonas NA",]
dedup <- dedup[rownames(dedup) != "Lysobacter tolerans",]

# Find the incatives
# inactive_tofind <- rownames(dedup)[rowSums(dedup) < 0.035]
# inactive_tofind
# inactive_df <- orgkey[grep(c(paste0(c("Kytococcus", word(inactive_tofind, 1, sep = " ")), collapse = "|")), orgkey$genus),] %>%
#   dplyr::filter(!grepl("avermitilis", orgs))

pdf("output/substrate_comparisons/substrate_comparison_heatmap_unclustered_C6_only_log10_scale.pdf", width = 4, height = 8)
pheatmap(
  cluster_cols = F,
  cluster_rows = T,
  border_color = NA,
  mat = dedup, 
  breaks = mat_breaks,
  color = pal2,
  annotation_names_row = T, 
  fontsize_col = 10, 
  fontsize_row = 8, 
  annotation_names_col =  T)
dev.off()
dedup[,1]

pdf("output/substrate_comparisons/substrate_comparison_heatmap_unclustered_C6_only_transposed_log10_scale.pdf", width = 12, height = 4)
pheatmap(
  cluster_cols = T,
  border_color = NA,
  mat = t(dedup), 
 #  breaks = mat_breaks,
  color = inferno(80),
  annotation_names_row = T, 
  fontsize_col = 8, 
  fontsize_row = 10, 
  annotation_names_col =  T)
dev.off()

