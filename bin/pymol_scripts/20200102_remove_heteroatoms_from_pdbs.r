# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest","gbm","nnet","e1071","svmpath","lars",
               "glmnet","svmpath", "data.table", "bio3d")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the PDB files
pdb <- read.pdb("data/caver_models/4ku5.pdb")
b.inds <- atom.select(pdb, chain="B", value = TRUE)
pdb_nohet <- atom.select(b.inds, "notprotein", inverse=TRUE, value = TRUE)

# Find starting coordinates

s143 <- atom.select(pdb_nohet, elety = c("CA"), resno = 143, value = TRUE)
s143$xyz[1,]
# -13.915  19.155 -26.472

write.pdb(pdb_nohet, "data/caver_models/4KU5_chainB_no_heteroatoms.pdb")


### Also remove heteroatoms from 3ROW
# Read in the PDB files
pdb <- read.pdb("data/caver_models/3row.pdb")
b.inds <- atom.select(pdb, chain="B", value = TRUE)
pdb_nohet <- atom.select(b.inds, "notprotein", inverse=TRUE, value = TRUE)
write.pdb(pdb_nohet, "data/caver_models/3row_chainB_no_heteroatoms.pdb")



### Also remove heteroatoms from X. oryzae 
pdb <- read.pdb("data/caver_models/3fk5.pdb")
# b.inds <- atom.select(pdb, chain="B", value = TRUE)
pdb_nohet <- atom.select(pdb, "notprotein", inverse=TRUE, value = TRUE)
write.pdb(pdb_nohet, "data/caver_models/3fk5_chainB_no_heteroatoms.pdb")

### Also remove heteroatoms from 4KTI
pdb <- read.pdb("data/caver_models/4kti.pdb")
b.inds <- atom.select(pdb, chain="A", value = TRUE)
pdb_nohet <- atom.select(b.inds, "notprotein", inverse=TRUE, value = TRUE)
write.pdb(pdb_nohet, "data/caver_models/4kti_chainA_no_heteroatoms.pdb")
