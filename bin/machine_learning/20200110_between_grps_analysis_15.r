# Install packages
pacman::p_load(DECIPHER, tidyverse, seqinr, bgafun, RColorBrewer)

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

### Extract residues within 12 Angstroms
sqs <- readAAStringSet("data/alignments/73_OleA_JGI_unaligned.fasta")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})

sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
source("lib/extract_12angstrom_residues.R")

extract_84_list <- lapply(1:length(sqs), function(x) { extract_12angstrom_residues(query_fils[x]) }) # takes time to run
extract_84_df <- data.frame(matrix(unlist(extract_84_list), nrow = length(extract_84_list), byrow=T),
                            stringsAsFactors=FALSE)
colnames(extract_84_df)[1] <- "seq"

which_feats <- read_csv("data/machine_learning/84_residues_12_angstroms_4KU5_S143.csv")

comb_res <- which_feats %>% pull()
comb_res

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

supp1 <- read_csv("output/Supplemental_table_1_nsub_avg_activity.csv")
broad <- supp1$`NCBI Accession`[1:20]
narrow <- supp1$`NCBI Accession`[54:73]


brdnms <- dtf$column_label[grepl(paste0(c(broad, "4KU5"), collapse = "|"), dtf$column_label)]
nrwnms <- dtf$column_label[grepl(paste0(narrow, collapse = "|"), dtf$column_label)]
dtf$column_label[grepl(paste0(c(broad, "4KU5"), collapse = "|"), dtf$column_label)] <- paste0(brdnms, "_broad")
dtf$column_label[grepl(paste0(narrow, collapse = "|"), dtf$column_label)] <- paste0(nrwnms, "_narrow")
fullaaset <- AAStringSet(fullseq)
names(fullaaset) <- dtf$column_label
# writeXStringSet(fullaaset, "data/machine_learning/73_OleA_12angstrom_signature.faa")

# Write a trimmed alignment with only broad and narrow
brd <- fullaaset[grep("broad", names(fullaaset))]
nrw <- fullaaset[grep("narrow", names(fullaaset))]
comb <- AAStringSet(c(brd, nrw))

# Write to file
writeXStringSet(comb, "data/machine_learning/40_OleA_broad_narrow_sub_spec_12angstrom.faa")

# Read the alignment
j2tr <- readAAStringSet("data/machine_learning/40_OleA_broad_narrow_sub_spec_12angstrom.faa")
rdaln <- read.alignment("data/machine_learning/40_OleA_broad_narrow_sub_spec_12angstrom.faa", format = "fasta")
rdaln$seq <- toupper(rdaln$seq)

# Convert alignemnt into a binary presence-abscence matrix 
amino <- convert_aln_amino(rdaln)

# Define groups: "positive" and "negative"
grps <- rownames(amino)
grps[grep("_broad", grps)] <- "broad"
grps[grep("_narrow", grps)] <- "narrow"
grps <- as.factor(grps)


# Remove gaps
amino.gapless <- remove_gaps_groups(x = amino, z = grps) 

# Run correspondance analysis (CA)
ca.aap <- bga(t(amino.gapless + 1), grps)

# Identify important residues
top_res <- top_residues_2_groups(ca.aap)
names(top_res) <- gsub("X", "", names(top_res))

# Conserved amino acid profiles
profiles <- create_profile_strings(amino.gapless, grps)
keep <- profiles[, colnames(profiles) %in% names(top_res)]

featdf <- which_feats %>%
  rownames_to_column() %>%
  dplyr::mutate(ind = as.numeric(rowname))

keepdf <- data.frame(keep, stringsAsFactors = F) %>%
  t() %>%
  data.frame() %>%
  rownames_to_column(.) %>%
  dplyr::mutate(aa = substr(rowname, nchar(rowname), nchar(rowname))) %>%
  dplyr::mutate(ind = as.numeric(gsub("[[:alpha:]]", "", word(rowname, -1, sep = "X")))) %>%
  left_join(., featdf, by = "ind")
keepdf

# Channel residues
imp_chan_res <- c(261, 284, 292, 203, 172, 173)
imp_non_chan_res <- c(320, 321, 318, 323, 313, 240, 310)
channel_a <- c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353)
channel_b <- c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246)

chana_df <- keepdf %>%
  dplyr::filter(ang12.resi %in% channel_a) # T292 is important (in 20 of the highly active ones)
chana_df  
table(substr(rdaln$seq, 55, 55))

chanb_df <- keepdf %>%
  dplyr::filter(ang12.resi %in% channel_b) # residue L203
chanb_df
table(substr(rdaln$seq, 29, 29))

nonchan <- keepdf %>%
  dplyr::filter(ang12.resi %in% imp_non_chan_res)
nonchan


