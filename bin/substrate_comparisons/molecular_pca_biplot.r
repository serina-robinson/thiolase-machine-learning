#library(devtools)
#install_github("vqv/ggbiplot")# Install packages
pacman::p_load("tidyverse", "ChemmineR", "readxl", "webchem", "Hmisc",
               "Rcpi", "recipes", "ggthemes", "caret", "earth", "ggbiplot")
# BiocManager::install("Rcpi", dependencies = c("Imports", "Enhances"))

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the raw data
rawdat <- read_csv("data/selected_molecular_properties_17pNPs.csv")

# Correlate with all activity 
activity <- read_csv("output/substrate_comparisons/enzyme_activity_per_substrate.csv")

# Combine the data with activity
mergdat <- rawdat %>% 
  na.omit() %>%
  left_join(., activity, by = c("cmpnd_abbrev" = "abbrev"))
rownames(mergdat) <- mergdat$cmpnd_abbrev

# Trim irrelevant features
trimdat <- mergdat %>%
  dplyr::select(-status, -cmpnd_abbrev, -IUPAC, -SMILES) # %>%
# as.data.frame()

# Split into training and testing
training_dat <- trimdat[!is.na(trimdat$activity),]
dim(training_dat)

testing_dat <- trimdat[is.na(trimdat$activity),]

# Remove zero-variance columns
novar_dat <- training_dat[ - as.numeric(which(apply(training_dat, 2, var) == 0))]

# Principal components analysis in R
pca_dat <- prcomp(novar_dat[,1:(ncol(novar_dat)-1)], center = TRUE, scale = TRUE)
pdf("output/PCA_biplot_PC1_PC2.pdf", width = 12, height = 12)
ggbiplot(pca_dat, obs.scale = 1, var.scale = 1) +
  scale_color_discrete(name = '') +
  theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

pca_dat$x[,1]

pdf("output/PCA_biplot_PC1_PC5.pdf", width = 12, height = 12)
ggbiplot(pca_dat, choices = c(1, 5), obs.scale = 1, var.scale = 1) +
  scale_color_discrete(name = '') +
  theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()
