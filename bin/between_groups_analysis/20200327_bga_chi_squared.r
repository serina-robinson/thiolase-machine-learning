# Install packages
pacman::p_load(DECIPHER, tidyverse, seqinr, bgafun, RColorBrewer)

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

### Group size 35 ####

which_feats <- read_csv("data/machine_learning/84_residues_12_angstroms_4KU5_S143.csv")

j2tr <- readAAStringSet("data/machine_learning/70_OleA_broad_narrow_sub_spec_12angstrom.faa")
length(j2tr)
rdaln <- read.alignment("data/machine_learning/70_OleA_broad_narrow_sub_spec_12angstrom.faa", format = "fasta")
rdaln$seq <- toupper(rdaln$seq)

# Convert alignemnt into a binary presence-abscence matrix 
amino <- convert_aln_amino(rdaln)
amino[1,]

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
  dplyr::filter(ang12.resi %in% channel_a) %>%
  dplyr::mutate(channel = "A")# T292 is important (in 20 of the highly active ones)
chana_df  # T292

table(substr(rdaln$seq, 55, 55))

chanb_df <- keepdf %>%
  dplyr::filter(ang12.resi %in% channel_b) %>%
  dplyr::mutate(channel = "B") # residue L203
chanb_df

chanall_df <- chana_df %>%
  bind_rows(chanb_df)
chanall_df

chantrim <- chanall_df %>%
  dplyr::filter(ang12.resi == 203 | ang12.resi == 292)
chantrim

# Fisher's exact test is required because of low sample sizes

# T292
chisq.test(T292_tbl)

T292_tbl <- as.table(matrix(c(27,8,13,22), nrow = 2))
colnames(T292_tbl) <- c("broad", "narrow")
rownames(T292_tbl) <- c("T292", "other")
T292_tbl
fisher.test(T292_tbl) # 0.001485

# L203
chisq.test(L203_tbl)

L203_tbl <- as.table(matrix(c(31, 4, 17, 18), nrow = 2))
colnames(L203_tbl) <- c("broad", "narrow")
rownames(L203_tbl) <- c("T292", "other")
fisher.test(L203_tbl)
# p-value = 0.0006
