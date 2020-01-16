# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", "randomcoloR", "dendextend", "phylogram", "patchwork",
               "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales", "phytools", "ggtree")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "calculated_slopes"
fils

ll <- fils[grepl(suffix, fils)]
ll <- ll[!grepl("BocPhe|furf|scratch|only|averaged|benzoate|round|2019-11-15_7Ph", ll)] # rep1|rep2|reps|round|
ll

rawdat <- tibble(filename = ll) %>%
  mutate(file_contents = purrr::map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(cmpnd =  paste0(#word(word(filename, 2, sep = "output\\/\\/"), 1, sep = "\\/"), "_",
                                stringr::word(stringr::word(filename, 2, sep = "_"), 1, sep = "_"))) %>%
  dplyr::mutate(cmpnd = case_when(grepl("C6", filename) ~ "hexanoate",
                                  grepl("7Ph", cmpnd) ~ "7Ph heptanoate",
                                  grepl("Azido", cmpnd) ~ "azido",
                                  grepl("Butoxy", cmpnd) ~ "butoxy",
                                  grepl("Biotin", cmpnd) ~ "biotin",
                                  grepl("decanoate rep-1", cmpnd) ~ "decanoate",
                                  grepl("decanoate rep-2", cmpnd) ~ "decanoate",
                                  grepl("decanoate ", cmpnd) ~ "decanoate",
                                  grepl("ClPh", cmpnd) ~ "ClPh propionate",
                                  TRUE ~ cmpnd)) 

# Pull the already averaged hexanoate controls
avged_ctrls <- rawdat %>%
  dplyr::filter(grepl("hexanoate", cmpnd))

newdat <- rawdat %>%
  dplyr::filter(!grepl("hexanoate", cmpnd))

maprdat_log <- newdat %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope)) 
maprdat_log

# Combine with c6
maprdat_long <- maprdat_log %>%
  bind_rows(., avged_ctrls) %>%
  dplyr::select(cmpnd, org, log_slope)

# Find the average across replicates 
maprdat_avg <- maprdat_long %>%
  dplyr::group_by(cmpnd, org) %>% 
  dplyr::summarise_each(funs(mean), log_slope)
maprdat_avg

maprdat_merg <- as.data.frame(maprdat_avg, stringsAsFactors = F)

# Convert to wide format
maprdat_wide <- reshape2::dcast(maprdat_merg, org ~ cmpnd, value.var = "log_slope") 

maprdat_wide[is.na(maprdat_wide)] <- 0

rawdat_mat <- maprdat_wide %>%
  dplyr::select(-org) %>%
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


# Remove duplicates and exceptions
dedup <- maprdat_mat[!duplicated(rownames(maprdat_mat)),]
dedup_sort <- dedup[order(rowSums(dedup), decreasing = T),]

# Add in benzoate
testr <- cbind(dedup_sort, rep(0, nrow(dedup_sort)))
colnames(testr)[ncol(testr)] <- "benzoate"

testdf <- data.frame(testr, stringsAsFactors = F)
testdf$org <- rownames(testdf)
testdf_long <- gather(testdf, key = substrate, value = activity, X7Ph.heptanoate:benzoate, factor_key = TRUE)
testdf_long$substrate <- gsub("X7", "7", testdf_long$substrate)
head(testdf_long)

# write_csv(testdf_long, "output/substrate_comparisons/20200109_all_cmpnds_avg_log_slopes_for_modeling.csv")

# First convert to a hclust object
tr <- "data/73_OleA_aligned_trimmed_fasttree_phylo.nwk"
fastree <- treeio::read.tree(tr)
ggtree_phylo <- treeio::as.phylo(fastree)
ggt <- ggtree(ggtree_phylo, branch.length = F)

pdf("output/ggtree_results.pdf")
ggt
dev.off()

# Read in the taxonomy
rawtax <- read_excel("data/OleA_taxonomic_classification.xlsx") 
xantho <- rawtax[grep("Xanthomonas", rawtax$organism),]
xantho$organism <- c("Xanthomonas_campestris_OleA")
xantho$genus <- "Xanthomonas_campestris"

tax <- rawtax  %>%
  bind_rows(xantho)
tail(tax)
tax$organism <- gsub(" ", "_", tax$organism)

ptax <- ggtree_phylo$tip.label
ptax

for(i in 1:length(tax$organism)) {
  ind <- grep(tax$organism[i], ptax)
  tax$acc[i] <- word(ptax[ind], sep = "\\.1", 1)
  ptax[ind] <- paste0(ptax[ind], "_", tax$class[i], "_", tax$phylum[i])
}

tax$acc[is.na(tax$acc)] <- "NP_635607"

# Reorder the data frame to match
dt2 <- data.frame(ptax, word(ptax, sep = "\\.1", 1))
colnames(dt2) <- c("label", "acc")

# Merge with the tax df
merg <- dt2 %>%
  inner_join(tax, by = "acc") %>%
  dplyr::mutate(labs = paste0(word(organism, sep = "_", 1), " ",
                word(organism, sep = "_", 2)))

# Fix discrepancies in merging names
pseudo1 <- rownames(maprdat_mat)[grep("Pseudoxanthomonas", rownames(maprdat_mat))]
pseudo2 <- merg$labs[grep("Pseudoxanthomonas", merg$labs)][1]
merg$labs[grep("Pseudoxanthomonas", merg$labs)] <- pseudo1

leif1 <- rownames(maprdat_mat)[grep("Leifsonia", rownames(maprdat_mat))]
leif2 <- merg$labs[grep("Leifsonia", merg$labs)][1]
merg$labs[grep("Leifsonia", merg$labs)] <- leif1

rownames(maprdat_mat)[!rownames(maprdat_mat) %in% merg$labs] # everything matches

ggtree_phylo$tip.label <- merg$labs
ult <- phytools::force.ultrametric(ggtree_phylo, method=c("nnls","extend"))
root_lab <- ult$tip.label[grep("Myco", ult$tip.label)]
rtd_tree <- root(ult, root_lab)

hcl2 <- ape::as.hclust.phylo(rtd_tree)



resmat <- testr[match(ggtree_phylo$tip.label, rownames(testr)),]


newnames <- lapply(
  rownames(resmat),
  function(x) bquote(italic(.(x))))


respht <- pheatmap(
  cluster_cols = T,
  cluster_rows = F,
  border_color = NA,
  mat = resmat,
  color = pal2,
  annotation_names_row = T, 
  # clustering_callback = cl_cb,
  fontsize_col = 10, 
  fontsize_row = 8, 
  annotation_names_col =  T,
  labels_row = as.expression(newnames))


pdf("output/substrate_comparisons/20200109_substrate_comparison_heatmap_unclustered_log10_scale_per_hr_reps_avged_phylo.pdf", width = 6, height = 11)
respht
dev.off()


# Now add the phylogenetic tree dendrogram
wts <- as.numeric(match(ggtree_phylo$tip.label, rownames(testr)))
wts

cl_cb <- function(hcl, mat){
  
  wts <- match(ggtree_phylo$tip.label, rownames(mat))
  hclust_olo <- reorder(hcl, wts)
  return(hclust_olo)
}

wts <- match(ggtree_phylo$tip.label, rownames(testr))
wts
reorder(respht_clust$tree_row, wts, agglo.fun = "max")


respht_clust <- pheatmap(
  cluster_cols = T,
  cluster_rows = T,
  border_color = NA,
  mat = testr,
  color = pal2,
  annotation_names_row = T, 
  clustering_callback = cl_cb,
  fontsize_col = 10, 
  fontsize_row = 8, 
  annotation_names_col =  T,
  labels_row = as.expression(newnames))

#rn <- factor(rownames(testr), levels = rownames(testr))
#rn
#gt <- factor(ggtree_phylo$tip.label, levels = ggtree_phylo$tip.label)
# rn[order(rn, gt)]


ord_cols <- respht$tree_col$labels[respht$tree_col$order]
colnames(testr)

test_ord <- testr[,match(ord_cols, colnames(testr))]
colnames(test_ord)

ggt <- ggtree(ggtree_phylo) + 
  geom_tiplab(aes(fontface = "italic"),
              align = TRUE, linetype = 2, offset = 0.04) + 
  xlim(NA, 25)

#pdf("output/ggtree_with_heatmap.pdf", height = 20, width = 50)


# Try ggtree and barplot
avg_act <- test_ord %>%
  rowMeans(.) %>%
  data.frame(., stringsAsFactors = F) %>%
  rownames_to_column(., var = "org")

colnames(avg_act) <- c("org", "avg_activity")

# Merge with taxonomy
mergtax <- avg_act %>%
  left_join(., merg, by = c("org" = "labs"))
head(mergtax)


# Read in the custom palette
clustkey <- read_csv("data/OleA_palette_key.csv")
clustkey$levs[clustkey$levs == "Green non-sulfur bacteria"] <- "Chloroflexi"

clustkey
levels(as.factor(ggdf$data$class))


mergtax$class[grepl("Opit", mergtax$class)] <- "Opitutae"
mergord <- mergtax[rev(match(ggt$data$label[ggt$data$isTip], mergtax$org)),] %>%
  dplyr::select(-label) %>%
  dplyr::mutate(label = org) %>%
  dplyr::left_join(., clustkey, by = c("class" = "levs")) %>%
  dplyr::select(label, org, avg_activity, acc, organism, genus, family, order, class, phylum, pal2)

mergord[is.na(mergord$pal2),]

ggt$panel <- 'Tree'
mergord$panel <- 'Stats'
ggt$panel <- factor(ggt$panel, levels=c("Tree", "Stats"))
mergord$panel <- factor(mergord$panel, levels=c("Tree", "Stats"))


pdf("output/heatmap_combo_barplot.pdf", width = 30, height = 20)
p <- ggplot(mapping=aes(x = org, 
                                        y = avg_activity,
                                        color = class,
                                        fill = class)) +
  facet_grid(.~panel, scale="free_y") + theme_tree2()
g2 <- p + geom_bar(data = mergord, stat= "identity") + 
  coord_flip() + 
  ggt
gheatmap(g2, test_ord, offset = 6.5, width = 2.5, font.size = 1, color = NULL) +
scale_fill_viridis(option = "inferno")
dev.off()



ggplot(data = mergord, mapping=aes(x = org, 
                                   y = avg_activity,
                                   color = class,
                                   fill = class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggt +
  facet_grid(.~panel, scale="free_x")



pdf("output/ggtree_with_barplot.pdf",width = 10, height = 20)
ggplot(data = avg_act, mapping=aes(x=x, y=y)) + facet_grid(.~panel, scale="free_x")
dev.off()

mergord$label


ggdf <- ggt %<+% mergord
ggdf$data$class

clustord <- clustkey[match(levels(as.factor(ggdf$data$class)), clustkey$levs),]

pdf("output/ggtree_with_circsize.pdf",width = 10, height = 20)
gsz <- ggdf +
  geom_tippoint(aes(size = avg_activity, color = class, fill = class, alpha = avg_activity), x = 20) +
  scale_fill_manual(values = clustord$pal2) +
  scale_color_manual(values = clustord$pal2)
gsz
dev.off()

ggdf$data$pal2
gsz <- ggdf +
  geom_tippoint(aes(size = avg_activity, alpha = avg_activity), color = ggdf$data$pal2[ggdf$data$isTip], 
                fill = ggdf$data$pal2[ggdf$data$isTip], x = 12.5)
gsz  
# scale_fill_manual(values = clustord$pal2) +
  # scale_color_manual(values = clustord$pal2)



pdf("output/ggtree_heatmap_combo_no_legend.pdf", width = 11, height = 13.75)
gheatmap(gsz, test_ord, offset = 8, width = 2.4, font.size = 1, color = NULL) +
  scale_fill_viridis(option = "inferno") +
  theme(legend.title = element_blank(), legend.position = "none")
dev.off()

# gheatmap(ggt, testr, low = pal2[1], high = pal2[70])


pdf("output/substrate_comparisons/20200109_substrate_comparison_heatmap_unclustered_transposed_log10_scale_per_hr_reps_avged.pdf", width = 15, height = 8.5)
pheatmap(
  cluster_cols = T,
  cluster_rows = T,
  border_color = NA,
 #  mat = t(dedup_sort), 
  mat = t(testr),
  #  breaks = mat_breaks,
  color = pal2,
  angle_col = 315,
  annotation_names_row = T, 
  fontsize_col = 9, 
  fontsize_row = 10, 
  annotation_names_col = T,
  labels_col = as.expression(newnames))
dev.off()



