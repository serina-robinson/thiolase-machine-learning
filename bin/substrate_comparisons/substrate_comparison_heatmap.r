# Pipeline for files with # Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales", "wesanderson")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "_calculated_slopes.csv"


ll <- fils[grepl(suffix, fils)]
rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(cmpnd = word(word(filename, 2, sep = "_"), 1, sep = "_"))

head(rawdat)

# Write to file
maprdat_long  <- rawdat
write_csv(maprdat_long, "output/substrate_comparisons/all_slopes_long_format.csv")



maprdat_wide <- dcast(maprdat_long, org ~ cmpnd, value.var = "max_slope") 
maprdat_wide[is.na(maprdat_wide)] <- 0

maprdat_mat <- maprdat_wide %>%
  dplyr::select(-org) %>%
  as.matrix()

summary(maprdat_long$max_slope)
        
maprdat_mat[is.na(maprdat_mat)]


rownames(maprdat_mat) <- maprdat_wide$org
maprdat_mat <- maprdat_mat #%>%
#  t()

quantile_breaks <- function(xs, n = 20) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(maprdat_mat, n = 21)
mat_breaks
mat_breaks <- sort(c(mat_breaks, 0.3, 0.5, 0.7))
mat_breaks
# labs <- as.character(scientific(mat_breaks, digit = 3))
maprdat_mat
pal <- inferno(length(mat_breaks))
pal2 <- pal[c(1, 3:11)]
pal2

tax <- read_excel("data/OleA_taxonomic_classification.xlsx")
rownames(maprdat_mat) <- gsub("_", " ", rownames(maprdat_mat))
rownames(maprdat_mat) <- gsub("XC", "Xanthomonas campestris OleA*", rownames(maprdat_mat))
rownames(maprdat_mat) <- gsub("Pseudoxanthomonas", "Pseudoxanthomonas sp.", rownames(maprdat_mat))
rownames(maprdat_mat) <- paste0(word(rownames(maprdat_mat), 1, sep = " "), " ", word(rownames(maprdat_mat), 2, sep = " "))

annot_pal <- colorRampPalette(brewer.pal(8, "Set2"))(8)
annot_pal <- c("dodgerblue", annot_pal[4:6])

# pdf("output/annot_pal.pdf")
# show_col(annot_pal)
# dev.off()

annot_df <- data.frame(maprdat_mat) %>%
  dplyr::mutate(orgs = rownames(.)) %>%
  dplyr::mutate(genus = word(orgs, 1, sep = " ")) %>%
  left_join(., tax, by = "genus") %>%
  dplyr::filter(!duplicated(orgs)) %>%
  dplyr::select(orgs, genus, phylum) %>%
  dplyr::mutate(fill = annot_pal[as.numeric(as.factor(phylum))])
dim(annot_df)

annot_nams <- data.frame(phylum = as.factor(annot_df$phylum))
annot_nams
annot_cols <- list(phylum = #c(levels(annot_nams$ID)[1] = annot_pal[1],
                   #levels(annot_nams$ID)[2] = annot_pal[2]))
                   c("Actinobacteria" = annot_pal[1],
                          "Firmicutes" = annot_pal[2],
                           "Proteobacteria" = annot_pal[3],
                     "Chlamydiae" = annot_pal[4]))
                            
annot_cols
rownames(maprdat_mat)
annot_nams

rownames(annot_nams) <- rownames(maprdat_mat)
newnames <- lapply(
  rownames(maprdat_mat),
  function(x) bquote(italic(.(x))))


pdf("output/substrate_comparisons/substrate_comparison_heatmap_clustered_only_active_no_labs.pdf", width = 9, height = 6)
pheatmap(annotation_row = annot_nams,
         annotation_colors = annot_cols[1],
         labels_col = rep("", 6),
         labels_row = as.expression(newnames),
         border_color = NA,
         mat = maprdat_mat, 
         breaks = mat_breaks,
         color = pal2,
         annotation_names_row = F, 
         fontsize_col = 10, 
         fontsize_row = 10, 
         annotation_legend = T)
dev.off()

# Ggplot2 heatmap
# pdf("output/substrate_comparison_heatmap.pdf")
# ggplot(maprdat_long, aes(cmpnd, variable)) +
#   geom_tile(aes(fill = value)) + 
#  #  geom_text(aes(label = round(value, 1))) +
#   scale_fill_gradient(low = "white", high = "royalblue") +
#   theme(axis.text.x = element_text(angle = 90))
# dev.off()

# Now combine with all organisms not present 
orgkey <- read_csv("data/72_OleA_masterwell_org_key.csv")
inactives <- orgkey$orgs[!orgkey$orgs %in% rownames(maprdat_mat)]
inactive_mat <- matrix(ncol = ncol(maprdat_mat), nrow = length(inactives), 0) #colnames = colnames(maprdat_mat))
inactive_mat
colnames(inactive_mat) <- colnames(maprdat_mat)
rownames(inactive_mat) <- inactives
inactive_mat

full_mat <- rbind(maprdat_mat, inactive_mat)
rownames(full_mat) <- gsub("_", " ", rownames(full_mat))
rownames(full_mat) <- paste0(word(rownames(full_mat), 1, sep = " "), " ", word(rownames(full_mat), 2, sep = " "))
rownames(full_mat) <- gsub("XC", "Xanthomonas campestris OleA*", rownames(full_mat))

pdf("output/substrate_comparisons/substrate_comparison_heatmap_clustered_full.pdf", width = 6, height = 10)
pheatmap(
  border_color = NA,
  full_mat, 
  breaks = mat_breaks,
  color = inferno(10),
  annotation_names_row = T, 
  fontsize_col = 10, 
  fontsize_row = 8, 
  annotation_names_col =  T)
dev.off()

# Choose active OleAs 
# Write interesting ones to file
gammas <- c("Xanthomonas", "Wenz", "Chromato", "Silani", "Areni", "Lutei", "Thermomon", "Pseudoxanthomon", 
            "XC", "Kyto", "Mobili", "Actinoplan", "Granulosic")
maprdat_active <- maprdat_long %>%
  dplyr::add_count(org) %>%
  dplyr::filter(n >= 2) %>%
  dplyr::filter(!grepl(paste0(gammas, collapse = "|"), org)) %>%
  inner_join(., orgkey, by = c("org" = "orgs")) %>%
  arrange(desc(max_slope))
unique(maprdat_active$org)
writeLines(word(unique(maprdat_active$gene), sep = " ", 1), "output/active_distant_OleA_accessions.txt")

maprdat_active
write_csv(maprdat_active, "output/substrate_comparisons/active_non_Xantho_hits.csv")

