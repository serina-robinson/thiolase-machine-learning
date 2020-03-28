# Install packages
pacman::p_load(DECIPHER, tidyverse, seqinr, bgafun, RColorBrewer)

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

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
# summary(supp1$`Average activity`)
# supp1$`Average activity`[20]
# supp1$`Average activity` > 0.8264
# broad <- supp1$`NCBI Accession`[1:20]
# narrow <- supp1$`NCBI Accession`[54:73]

# Write Group size 25 to file
broad <- supp1$`NCBI Accession`[1:25]
narrow <- supp1$`NCBI Accession`[49:73]

brdnms <- dtf$column_label[grepl(paste0(c(broad, "4KU5"), collapse = "|"), dtf$column_label)]
nrwnms <- dtf$column_label[grepl(paste0(narrow, collapse = "|"), dtf$column_label)]
dtf$column_label[grepl(paste0(c(broad, "4KU5"), collapse = "|"), dtf$column_label)] <- paste0(brdnms, "_broad")
dtf$column_label[grepl(paste0(narrow, collapse = "|"), dtf$column_label)] <- paste0(nrwnms, "_narrow")
fullaaset <- AAStringSet(fullseq)
names(fullaaset) <- dtf$column_label
names(fullaaset)
# writeXStringSet(fullaaset, "data/machine_learning/73_OleA_12angstrom_signature.faa")

# Write a trimmed alignment with only broad and narrow
broad[broad == "NP_635607"] <- "4KU5"
brd <- fullaaset[grep(paste0(broad, collapse = "|"), names(fullaaset))]

nrw <- fullaaset[grep(paste0(narrow, collapse = "|"), names(fullaaset))]
length(nrw)

comb <- AAStringSet(c(brd, nrw))
length(comb)

# Write Group size 25 to file
writeXStringSet(comb, "data/machine_learning/50_OleA_broad_narrow_sub_spec_12angstrom.faa")


###  Group size 35
broad <- supp1$`NCBI Accession`[1:35]
narrow <- supp1$`NCBI Accession`[39:73]
length(narrow)

brdnms <- dtf$column_label[grepl(paste0(c(broad, "4KU5"), collapse = "|"), dtf$column_label)]
nrwnms <- dtf$column_label[grepl(paste0(narrow, collapse = "|"), dtf$column_label)]
dtf$column_label[grepl(paste0(c(broad, "4KU5"), collapse = "|"), dtf$column_label)] <- paste0(brdnms, "_broad")
dtf$column_label[grepl(paste0(narrow, collapse = "|"), dtf$column_label)] <- paste0(nrwnms, "_narrow")
fullaaset <- AAStringSet(fullseq)
names(fullaaset) <- dtf$column_label
names(fullaaset)
# writeXStringSet(fullaaset, "data/machine_learning/73_OleA_12angstrom_signature.faa")

# Write a trimmed alignment with only broad and narrow
broad[broad == "NP_635607"] <- "4KU5"
brd <- fullaaset[grep(paste0(broad, collapse = "|"), names(fullaaset))]

nrw <- fullaaset[grep(paste0(narrow, collapse = "|"), names(fullaaset))]
length(nrw)

comb <- AAStringSet(c(brd, nrw))
length(comb)

# Write Group size 15 to file
writeXStringSet(comb, "data/machine_learning/70_OleA_broad_narrow_sub_spec_12angstrom.faa")

## Group size 30
broad <- supp1$`NCBI Accession`[1:30]
narrow <- supp1$`NCBI Accession`[44:73]

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

# Write Group size 30 to file
# writeXStringSet(comb, "data/machine_learning/60_OleA_broad_narrow_sub_spec_12angstrom.faa")



# Write Group size 25 to file
broad <- supp1$`NCBI Accession`[1:25]
narrow <- supp1$`NCBI Accession`[49:73]

brdnms <- dtf$column_label[grepl(paste0(c(broad, "4KU5"), collapse = "|"), dtf$column_label)]
nrwnms <- dtf$column_label[grepl(paste0(narrow, collapse = "|"), dtf$column_label)]
dtf$column_label[grepl(paste0(c(broad, "4KU5"), collapse = "|"), dtf$column_label)] <- paste0(brdnms, "_broad")
dtf$column_label[grepl(paste0(narrow, collapse = "|"), dtf$column_label)] <- paste0(nrwnms, "_narrow")
fullaaset <- AAStringSet(fullseq)
names(fullaaset) <- dtf$column_label
names(fullaaset)
# writeXStringSet(fullaaset, "data/machine_learning/73_OleA_12angstrom_signature.faa")

# Write a trimmed alignment with only broad and narrow
broad[broad == "NP_635607"] <- "4KU5"
brd <- fullaaset[grep(paste0(broad, collapse = "|"), names(fullaaset))]

nrw <- fullaaset[grep(paste0(narrow, collapse = "|"), names(fullaaset))]
length(nrw)

comb <- AAStringSet(c(brd, nrw))
length(comb)

# Write Group size 25 to file
writeXStringSet(comb, "data/machine_learning/50_OleA_broad_narrow_sub_spec_12angstrom.faa")


###  Group size 20 is done

### Group size 15
broad <- supp1$`NCBI Accession`[1:15]
narrow <- supp1$`NCBI Accession`[59:73]
length(broad)

brdnms <- dtf$column_label[grepl(paste0(c(broad, "4KU5"), collapse = "|"), dtf$column_label)]
nrwnms <- dtf$column_label[grepl(paste0(narrow, collapse = "|"), dtf$column_label)]
dtf$column_label[grepl(paste0(c(broad, "4KU5"), collapse = "|"), dtf$column_label)] <- paste0(brdnms, "_broad")
dtf$column_label[grepl(paste0(narrow, collapse = "|"), dtf$column_label)] <- paste0(nrwnms, "_narrow")
fullaaset <- AAStringSet(fullseq)
names(fullaaset) <- dtf$column_label
names(fullaaset)
# writeXStringSet(fullaaset, "data/machine_learning/73_OleA_12angstrom_signature.faa")

# Write a trimmed alignment with only broad and narrow
broad[broad == "NP_635607"] <- "4KU5"
brd <- fullaaset[grep(paste0(broad, collapse = "|"), names(fullaaset))]

nrw <- fullaaset[grep(paste0(narrow, collapse = "|"), names(fullaaset))]
length(nrw)

comb <- AAStringSet(c(brd, nrw))
length(comb)

# Write Group size 15 to file
writeXStringSet(comb, "data/machine_learning/30_OleA_broad_narrow_sub_spec_12angstrom.faa")

