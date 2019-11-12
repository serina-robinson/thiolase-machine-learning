# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", "randomcoloR", "RColorBrewer", 
               "ggplot2", "pheatmap", "viridis", "scales")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "calculated_slopes"
fils

ll0 <- fils[grepl(suffix, fils)]
ll0
ll <- ll0[grepl("hexanoate|C6", ll0)]
ll
rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(date_run = substr(word(filename, 1, sep = "_"), nchar("output//2019-11-08/") + 1, nchar("output//2019-11-08/") + 10)) %>%
  dplyr::mutate(cmpnd = case_when(grepl("round1", filename) ~ "hexanoate round 1",
                                  grepl("round2", filename) ~ "hexanoate round 2",
                                  TRUE ~ "hexanoate round 3"))

table(rawdat$cmpnd)

# Write to file
maprdat_long  <- rawdat 
dat_order <- unique(rawdat$cmpnd[order(rawdat$date_run)])
dat_order
maprdat_merg <- as.data.frame(maprdat_long, stringsAsFactors = F)

maprdat_wide <- reshape2::dcast(maprdat_merg, org ~ cmpnd, value.var = "max_slope") 
maprdat_wide[is.na(maprdat_wide)] <- 0

# Order by date


maprdat_mat <- maprdat_wide %>%
  dplyr::select(-org) %>%
  dplyr::select(dat_order) %>%
  as.matrix()


rownames(maprdat_mat) <- maprdat_wide$org
maprdat_mat <- maprdat_mat #%>%
colnames(maprdat_mat)
#  t()

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
rownames(maprdat_mat)

annot_pal <- colorRampPalette(brewer.pal(8, "Set2"))(8)
annot_pal <- c("dodgerblue", annot_pal[4:6], "purple", "red", "orange")
annot_pal

annot_df <- data.frame(maprdat_mat) %>%
  dplyr::mutate(orgs = rownames(.)) %>%
  dplyr::mutate(genus = word(orgs, 1, sep = " ")) %>%
  left_join(., tax, by = "genus") %>%
  dplyr::filter(!duplicated(orgs)) %>%
  dplyr::select(orgs, genus, phylum) %>%
  dplyr::mutate(fill = annot_pal[as.numeric(as.factor(phylum))])
dim(annot_df)

annot_df$phylum[grep("Verruco", annot_df$phylum)] <- "Verrucomicrobia"
annot_df$phylum
annot_nams <- data.frame(phylum = as.factor(annot_df$phylum))

annot_cols <- list(phylum = # c(levels(annot_nams$phylum)[1] = annot_pal[1],
                              # levels(annot_nams$phylum)[2] = annot_pal[2],
                              # levels(annot_nams$phylum)[3] = annot_pal[3],
                              # levels(annot_nams$phylum)[4] = annot_pal[4],
                              # levels(annot_nams$phylum)[5] = annot_pal[5],
                              # levels(annot_nams$phylum)[6] = annot_pal[6]))
                   #c(levels(annot_nams$ID)[1] = annot_pal[1],
                   #levels(annot_nams$ID)[2] = annot_pal[2]))
                   c("Actinobacteria" = annot_pal[1],
                     "Firmicutes" = annot_pal[2],
                     "Proteobacteria" = annot_pal[3],
                     "Chlamydiae" = annot_pal[4],
                     "Verrucomicrobia" = annot_pal[5],
                     "Chloroflexi" = annot_pal[6],
                     "Bacteroidetes" = annot_pal[7]))
                            
newnames <- lapply(
  rownames(maprdat_mat),
  function(x) bquote(italic(.(x))))

maprdat_fix <- maprdat_mat
rownames(maprdat_fix) <- substr(rownames(maprdat_mat), 1, 3)
sub_act <- data.frame((colSums(maprdat_fix)))
sub_act$abbrev <- rownames(sub_act)
colnames(sub_act) <- c( "activity", "abbrev")

# write_csv(sub_act, "output/substrate_comparisons/enzyme_activity_per_substrate.csv")
# pdf("output/activity_rowsums.pdf", width = 45, height = 10)
# barplot(sort(rowSums(maprdat_fix), decreasing = T), 
#         xlab = "organism",
#         ylab = "sum of activity across all substrates")
# dev.off()
# 
# pdf("output/activity_colsums.pdf", width = 15, height = 5)
# barplot(sort(colSums(maprdat_fix), decreasing = T), 
#         xlab = "substrate",
#         ylab = "sum of activity across all enzymes")
# dev.off()

# pdf("output/substrate_comparisons/substrate_comparison_heatmap_clustered_only_active_no_labs.pdf", width = 10, height = 10)
# pheatmap(annotation_row = annot_nams,
#          annotation_colors = annot_cols,
#          # labels_col = rep("", 6),
#          labels_row = as.expression(newnames),
#          border_color = NA,
#          mat = maprdat_mat, 
#          breaks = mat_breaks,
#          color = pal2,
#          annotation_names_row = F, 
#          fontsize_col = 10, 
#          fontsize_row = 10, 
#          annotation_legend = T)
# dev.off()

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

full_mat <- rbind(maprdat_mat, inactive_mat)

rownames(full_mat) <- gsub("_", " ", rownames(full_mat))
rownames(full_mat) <- paste0(word(rownames(full_mat), 1, sep = " "), " ", word(rownames(full_mat), 2, sep = " "))
dedup <- full_mat[!duplicated(rownames(full_mat)),]

# rownames(dedup)[rowSums(dedup) < 0.03]
dedup_df <- data.frame(dedup, stringsAsFactors = F)

dedup <- dedup[rownames(dedup) != "Pseudoxanthomonas NA",]
dedup <- dedup[rownames(dedup) != "Lysobacter tolerans",]
inactive_tofind <- rownames(dedup)[rowSums(dedup) < 0.035]
inactive_tofind

inactive_df <- orgkey[grep(c(paste0(c("Kytococcus", word(inactive_tofind, 1, sep = " ")), collapse = "|")), orgkey$genus),] %>%
  dplyr::filter(!grepl("avermitilis", orgs))

pdf("output/substrate_comparisons/substrate_comparison_heatmap_unclustered_C6_only.pdf", width = 4, height = 8)
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

pdf("output/substrate_comparisons/substrate_comparison_heatmap_unclustered_C6_only_transposed.pdf", width = 12, height = 4)
pheatmap(
  cluster_cols = T,
  border_color = NA,
  mat = t(dedup), 
  breaks = mat_breaks,
  color = pal2,
  annotation_names_row = T, 
  fontsize_col = 8, 
  fontsize_row = 10, 
  annotation_names_col =  T)
dev.off()

# # Heatmap with cutoff
cutoff <- 0.1
mat_breaks2 <- mat_breaks
mat_breaks2
mat_breaks2[mat_breaks2 <= 0.05] <- 0
mat_breaks2 <- unique(mat_breaks2)
# mat_breaks2
# palr <- inferno(length(mat_breaks))
#palr2 <- pal[c(1, 3:length(pal))]
# mat_breaks <- a2s.numeric(gsub(0.0001, cutoff, mat_breaks))
#dedup <- mapply(dedup, 2, 
                #function(x) { gsub(x < 0.1, 0, x) })
#head(dedup)

pdf(paste0("output/substrate_comparisons/substrate_comparison_heatmap_clustered_", cutoff, "_cutoff.pdf"), width = 6, height = 9)
pheatmap(
  border_color = NA,
  mat = dedup, 
  breaks = mat_breaks2,
  color = pal2,
  annotation_names_row = T, 
  fontsize_col = 10, 
  fontsize_row = 8, 
  annotation_names_col =  T)
dev.off()



# Choose active OleAs 
# Write interesting ones to file
gammas <- c("Xanthomonas", "Wenz", "Chromato", "Silani", "Areni", "Lutei", "Thermomon", "Pseudoxanthomon", 
            "XC", "Kyto", "Mobili", "Actinoplan", "Granulosic")
maprdat_active <- maprdat_merg %>%
  dplyr::add_count(org) %>%
  dplyr::filter(n >= 2) %>%
  dplyr::filter(!grepl(paste0(gammas, collapse = "|"), org)) %>%
  inner_join(., orgkey, by = c("org" = "orgs")) %>%
  arrange(desc(max_slope))
unique(maprdat_active$org)
writeLines(word(unique(maprdat_active$gene), sep = " ", 1), "output/active_distant_OleA_accessions.txt")

maprdat_active
write_csv(maprdat_active, "output/substrate_comparisons/active_non_Xantho_hits.csv")

