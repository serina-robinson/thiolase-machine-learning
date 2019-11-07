# Install packages
pacman::p_load("tidyverse", "Biostrings", "DECIPHER")

# Read in the sequences
aa <- readAAStringSet("data/alignments/73_OleA_JGI_unaligned.fasta")
aa.al <- AlignSeqs(aa)
BrowseSeqs(aa.al)
aa.sub <- AAStringSet(substr(aa.al, start = 25, stop = 459))
matchPattern("NAS", aa.sub$`4KU5_Xanthomonas_campestris`)
aa.sub$`4KU5_Xanthomonas_campestris`[159] <- "C"
aa.sub$`4KU5_Xanthomonas_campestris`[159]

writeXStringSet(aa.sub, "data/alignments/73_OleA_aligned_trimmed_4KU5.fasta")
