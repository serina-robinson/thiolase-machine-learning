# Install packages
pacman::p_load("tidyverse", "ChemmineR", "readxl", "webchem", "Rcpi")
# BiocManager::install("Rcpi", dependencies = c("Imports", "Enhances"))
# BiocManager::install("ChemmineOB", dependencies = c("Imports", "Enhances"))

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the spreadsheet of molecules and SMILES
cmpnds <- read_excel("data/2019-11-01/cmpnds_tested_as_of_2019-11-01.xlsx")
# cids <- get_cid(cmpnds$IUPAC) # doesn't work because not in pNP
cids
smiles <- cmpnds$SMILES
smiles

# Read in the SDFset
sdfset <- read.SDFset("data/17_pNPs.sdf")

propma <- atomcountMA(sdfset[1:length(sdfset)], addH=TRUE) 
boxplot(propma, col="blue", main="Atom Frequency") 

# Calculate Chemmine properties
propdf <- data.frame(MF=MF(sdfset, addH=TRUE), MW=MW(sdfset, addH=FALSE),
                     Ncharges=sapply(bonds(sdfset, type="charge"), length),
                     atomcountMA(sdfset, addH=FALSE), 
                     groups(sdfset, type="countMA"), 
                     rings(sdfset, upper=6, type="count", arom=TRUE))
propdf
# Calculate Open Babel properties
sdffil <- "data/17_pNPs.sdf"
mol <- readMolFromSDF(sdffil)
# aliphat <- extractDrugLongestAliphaticChain(mol)
all_props <- extractDrugAIO(mol) # calculates all properties at once



# write_csv(full_df, "data/all_molecular_properties.csv")

dup_cols_bool <- apply(full_df, 2, function(x) all(duplicated(x)))
dup_cols_bool

all_prop_df <- data.frame(all_props) %>%
  bind_cols(propdf) %>%
  dplyr::select(-MF) %>%
  select_if(~ !all(is.na(.))) %>%
  select_if(~ !all(duplicated(.))) %>%
  select_if(colSums(.) != 0) # 142 of them are relevant

full_df <- cmpnds %>%
  dplyr::mutate(MW = propdf$MW) %>%
  bind_cols(all_prop_df)
     
# write_csv(full_df, "data/selected_molecular_properties_17pNPs.csv")
#writeLines(colnames(full_df), sep = ",", "output/162_molecular_properties.txt")
