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

testm <- testr
testm[testm > 0] <- 1
rownames(testm)
nsub <- data.frame(cbind(rownames(testm), rowSums(testm)))
colnames(nsub) <- c("org", "nsub")


# First convert to a hclust object
tr <- "data/73_OleA_aligned_trimmed_fasttree_phylo.nwk"
fastree <- treeio::read.tree(tr)
ggtree_phylo <- treeio::as.phylo(fastree)
ggt <- ggtree(ggtree_phylo, branch.length = F)

# Read in the taxonomy
rawtax <- read_excel("data/OleA_taxonomic_classification.xlsx") 
table(rawtax$phylum)
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

# Create a phylogenetic tree
ggtree_phylo$tip.label <- merg$labs
namord <- names(sort(colMeans(testr), decreasing = T))
namord
test_ord <- testr[,match(namord, colnames(testr))]

ggt <- ggtree(ggtree_phylo) + 
  geom_tiplab(aes(fontface = "italic"),
              align = TRUE, linetype = 2, offset = 0.04) + 
  xlim(NA, 20)
ggt

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

mergtax$class[grepl("Opit", mergtax$class)] <- "Opitutae"
mergord <- mergtax[rev(match(ggt$data$label[ggt$data$isTip], mergtax$org)),] %>%
  dplyr::select(-label) %>%
  dplyr::mutate(label = org) %>%
  dplyr::left_join(., clustkey, by = c("class" = "levs")) %>%
  dplyr::select(label, org, avg_activity, acc, organism, genus, family, order, class, phylum, pal2)


df2 <- data.frame(cbind(mergord$acc, mergord$organism, mergord$label, as.numeric(mergord$avg_activity)), stringsAsFactors = F)
df2

df3 <- df2[order(df2$X2, decreasing = T),] %>%
  dplyr::left_join(., nsub, by = c("X3" = "org")) %>%
  dplyr::mutate(organism = gsub("_", " ", X2)) %>%
  dplyr::mutate(activity = round(as.numeric(X4), 3)) %>%
  dplyr::select(X1, organism, activity, nsub) %>%
  arrange(desc(activity))

df3
df3[df3$`Average activity` > 1,] # 25 enzymes with activity greater than 1
nrow(df3[df3$`Average activity` < 0.55,]) # 25 enzymes with activity less than 0.55
df3[df3$`Average activity` < 0.55,]

colnames(df3) <- c('NCBI Accession', 'Organism', 'Average activity', 'Total number of substrates accepted')
# write_csv(df3, "output/Supplemental_table_1_nsub_avg_activity.csv")

ggdf <- ggt %<+% mergord

clustord <- clustkey[match(levels(as.factor(ggdf$data$class)), clustkey$levs),]

pdf("output/ggtree_with_circsize.pdf",width = 10, height = 20)
gsz <- ggdf +
  geom_tippoint(aes(size = avg_activity, color = class, fill = class), x = 20) +
  scale_fill_manual(values = clustord$pal2) +
  scale_color_manual(values = clustord$pal2)
gsz
dev.off()

gsz <- ggdf +
  geom_tippoint(aes(size = avg_activity), color = ggdf$data$pal2[ggdf$data$isTip], 
                fill = ggdf$data$pal2[ggdf$data$isTip], x = 10.35)
gsz  
# scale_fill_manual(values = clustord$pal2) +
  # scale_color_manual(values = clustord$pal2)



pdf("output/ggtree_heatmap_combo_no_legend.pdf", width = 11, height = 13.75)
gheatmap(gsz, test_ord, offset = 5.75, width = 2, font.size = 1, color = NULL) +
  scale_fill_viridis(option = "inferno") +
  theme(legend.title = element_blank(), legend.position = "none")
dev.off()

pdf("output/ggtree_heatmap_combo_with_legend.pdf", width = 11, height = 13.75)
gheatmap(gsz, test_ord, offset = 5.75, width = 1.75, font.size = 1, color = NULL) +
  scale_fill_viridis(option = "inferno") +
  theme(legend.title = element_blank())
dev.off()






