# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggseqlogo")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the reference
sqs <- readAAStringSet("data/alignments/73_OleA_JGI_unaligned.fasta")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})


sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
source("lib/extract_12angstrom_residues.R")

extract_84_list <- lapply(1:length(sqs), function(x) { extract_12angstrom_residues(query_fils[x]) })

extract_84_df <- data.frame(matrix(unlist(extract_84_list), nrow = length(extract_84_list), byrow=T), 
                            stringsAsFactors=FALSE)
colnames(extract_84_df)[1] <- "seq"



chan_res <- c(261, 284, 292, 203, 172, 173)
non_chan_res <- c(320, 321, 318, 323, 313, 240, 310)
comb_res <- c(chan_res, non_chan_res)
comb_res
which_feats <- read_csv("data/machine_learning/84_residues_12_angstroms_4KU5_S143.csv")

resll <- list()
for(i in 1:length(comb_res)) {
  ind <- grep(comb_res[i], which_feats$ang12.resi)
  resll[[i]] <- substr(extract_84_df$seq, ind, ind)
}

# 13 residues of interest
names(resll) <- paste0("X", comb_res)
dtf <- dplyr::bind_rows(resll, .id = "column_label")
dtf$column_label <- names(sqs)
dtf$fullseq <- unite(dtf, "fullseq", 2:ncol(dtf), sep = "")
fullseq <- pull(dtf$fullseq)
fullseq

pdf("output/important_res_logo.pdf", width = 4, height = 5)
p <- ggplot() + geom_logo(fullseq, method = "p", col_scheme = 'chemistry') + theme_logo() + 
  theme(legend.position = 'none', axis.text.x = element_blank()) 
p
dev.off()

just_super_active <- dtf %>%
  dplyr::filter(grepl("Kytococcus|Mobilicoccus|Granulosicoccus|Actinoplanes|Thermomonas|Halo", column_label))
active_fullseq <- pull(just_super_active$fullseq)

pdf("output/active_enzymes_res_logo.pdf", width = 4, height = 2.5)
p <- ggplot() + geom_logo(active_fullseq, method = "p", col_scheme = 'chemistry') + theme_logo() + 
  theme(legend.position = 'none', axis.text.x = element_blank()) 
p
dev.off()

xantho_only <- dtf %>%
  dplyr::filter(grepl("4KU5", column_label))
xantho_only
xantho_only <- pull(xantho_only$fullseq)
xantho_only

pdf("output/4KU5_important_res_logo.pdf", width = 4, height = 1.5)
p <- ggplot() + geom_logo(xantho_only, method = "p", col_scheme = 'chemistry') + theme_logo() + 
  theme(legend.position = 'none', axis.text.x = element_blank()) 
p
dev.off()
