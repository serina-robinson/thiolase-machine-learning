# Install packages
pacman::p_load("tidyverse", "ChemmineR", "readxl", "webchem", "Rcpi")
# BiocManager::install("Rcpi", dependencies = c("Imports", "Enhances"))
# BiocManager::install("ChemmineOB", dependencies = c("Imports", "Enhances"))

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the spreadsheet of molecules and SMILES
cmpnds <- read_excel("data/substrate_comparisons/final_15_cmpnds_tested.xlsx")
smiles <- cmpnds$SMILES
writeLines(smiles) # To convert to SDF using OpenBabel
# http://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html

# Read in the SDFset
sdfset <- read.SDFset("data/substrate_comparisons/15_pNPs_openbabel_sdf_2D.sdf")

propma <- atomcountMA(sdfset[1:length(sdfset)], addH=TRUE) 
boxplot(propma, col="blue", main="Atom Frequency") 

# Calculate Chemmine properties
propdf <- data.frame(MF=MF(sdfset, addH=TRUE), MW=MW(sdfset, addH=FALSE),
                     Ncharges=sapply(bonds(sdfset, type="charge"), length),
                     atomcountMA(sdfset, addH=FALSE), 
                     groups(sdfset, type="countMA"), 
                     rings(sdfset, upper=6, type="count", arom=TRUE))
rownames(propdf) <- cmpnds$cmpnd_abbrev
propdf

# Calculate Open Babel properties
sdffil <- "data/substrate_comparisons/15_pNPs_openbabel_sdf_2D.sdf"
mol <- readMolFromSDF(sdffil)
# all_props <- extractDrugAIO(mol) # calculates all properties at once
raw_df <- all_props

write_csv(raw_df, "data/substrate_comparisons/15pNPs_293_molecular_properties.csv")

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
     
write_csv(full_df, "data/substrate_comparisons/15pNPs_159_selected_molecular_properties.csv")
#writeLines(colnames(full_df), sep = ",", "output/162_molecular_properties.txt")
