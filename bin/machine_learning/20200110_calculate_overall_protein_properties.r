# Install packages
pacman::p_load(Peptides, protr, DECIPHER, viridis, ggthemes, treeio, tidyverse, seqinr, bgafun,
               RColorBrewer, ape, phangorn, Biostrings, data.table, purrr, dplyr, ggtree)

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

### Read in full-length sequences
sqs <- readAAStringSet("data/alignments/73_OleA_JGI_unaligned.fasta")

### Use wild-type OleA rather than 4KU5
olea <- readAAStringSet("data/machine_learning/Xc_OleA.faa")

# Combine and write to file
comb <- AAStringSet(c(olea, sqs[!grepl("4KU5", names(sqs))]))
writeXStringSet(comb, "data/alignments/73_OleA_with_wildtype_Xc.fasta")

# Read in full-length sequences
sqs <-read.alignment("data/alignments/73_OleA_with_wildtype_Xc.fasta", format = "fasta")
sqs$seq <- toupper(sqs$seq)

# Read in 12 Å core
coreaa <- readAAStringSet("data/machine_learning/73_OleA_12angstrom_signature.faa")
core <- read.alignment("data/machine_learning/73_OleA_12angstrom_signature.faa", format = "fasta")
core$seq <- toupper(core$seq)
  
## Using the Peptides package
# Cruciani 
crciani <- sapply(sqs$seq, crucianiProperties)
cdf <- do.call(rbind.data.frame, crciani)
colnames(cdf) <- names(crciani[[1]])

# Fasgai
fsgai <- sapply(sqs$seq, fasgaiVectors)
fdf <- do.call(rbind.data.frame, fsgai)
colnames(fdf) <- names(fsgai[[1]])
  
# ST scale
st_scale <- sapply(sqs$seq, stScales)
stdf <- do.call(rbind.data.frame, st_scale)
colnames(stdf) <- names(st_scale[[1]])

# T-scale
t_scale <- sapply(sqs$seq, tScales)
tdf <- do.call(rbind.data.frame, t_scale)
colnames(tdf) <- names(t_scale[[1]])

# VHSE
vhse_scale <- sapply(sqs$seq, vhseScales)
vdf <- do.call(rbind.data.frame, vhse_scale)
colnames(vdf) <- names(vhse_scale[[1]])

# Z-scale
z_scale <- sapply(sqs$seq, zScales)
zdf <- do.call(rbind.data.frame, z_scale)
colnames(zdf) <- names(z_scale[[1]])

# Kidera
kdera <- sapply(sqs$seq, kideraFactors)
kdf <- do.call(rbind.data.frame, kdera)
colnames(kdf) <- names(kdera[[1]])

# Combine everything into one data frame
sqdf <- data.frame(cbind(sqs$nam, sqs$seq, core$seq), stringsAsFactors = F)
colnames(sqdf) <- c("nams", "sqs", "core")
sqs$nam

combdf <- sqdf %>%
  dplyr::mutate(acc = word(nams, sep = "\\.1", 1)) %>%
  dplyr::mutate(aliphat_full = aIndex(sqs)) %>%
  dplyr::mutate(aliphat_core = aIndex(core)) %>%
  dplyr::mutate(hydrophob_full = hydrophobicity(sqs)) %>%
  dplyr::mutate(hydrophob_core = hydrophobicity(core)) %>%
  dplyr::mutate(boman_interactions = boman(sqs)) %>%
  dplyr::mutate(hmoment = hmoment(sqs)) %>%
  dplyr::mutate(instab = instaIndex(sqs)) %>%
  dplyr::mutate(lngth_prot = lengthpep(sqs)) %>%
  dplyr::mutate(mw_prot = mw(sqs)) %>%
  dplyr::mutate(pi_prot = pI(sqs)) %>%
  bind_cols(., zdf, tdf, cdf, vdf, stdf, fdf, kdf)
combdf

# write_csv(combdf, "data/machine_learning/73_overall_calculated_protein_properties.csv")



