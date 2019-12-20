# Install packages
pacman::p_load("tidyverse", "ChemmineR", "readxl", "webchem", "Rcpi", 
               "recipes", "ggthemes", "caret", "earth", "factoextra", "FactoMineR")
# BiocManager::install("Rcpi", dependencies = c("Imports", "Enhances"))

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the raw data
rawdat <- read_csv("data/selected_molecular_properties_17pNPs.csv") %>%
  dplyr::filter(!grepl("sulfinyl|furan", cmpnd_abbrev))
dim(rawdat)

trimdat <- rawdat %>%
  dplyr::select(-status, -cmpnd_abbrev, -IUPAC, -SMILES) # %>%
rownames(trimdat) <- rawdat$cmpnd_abbrev

# Only keep variables with nonzero variance
nozdat <- nearZeroVar(trimdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE]
which_rem

findat <- trimdat %>%
  dplyr::select(-which_rem)
head(findat)

# scaldat <- data.frame(scale(findat, center = T, scale = T), stringsAsFactors = F)
# head(scaldat)

# Principal components analysis in R
pca_dat <- prcomp(findat, center = TRUE, scale. = TRUE)
pca_abs <- abs(pca_dat$rotation)
head(pca_dat$rotation)
head(pca_abs)

res <- lapply(1:ncol(pca_dat$rotation), function(x) { rownames(pca_abs)[order(pca_abs[,x], decreasing = T)] })
names(res) <- paste0("PC", 1:15)

res_df <- bind_rows(res)
# write_csv(res_df, 'output/PC_loadings_molecular_descriptors_new.csv')
# for(i in 1:ncol(pca_dat$rotation)) {
#   res[[i]] <- names(pca_dat$rotation[order(pca_dat$rotation[,x], decreasing = T)]
# }
pdf("output/molecular_descriptor_PC_screeplot_variances.pdf")
screeplot(pca_dat, npcs = 15, type = "lines")
dev.off()

pca_dat$x
findat
dat <- data.frame(pca_dat$x[,1:2], stringsAsFactors = F)
dat


pdf(paste0("output/molecular_descriptor_pca_dat.pdf"), width = 5, height = 5)
par(mar=c(0.01,0.01,0.01,0.01))
ggplot(data = dat, aes(x = PC1, y = PC2)) + 
  geom_text(aes(label = rownames(dat))) +
  scale_shape(solid = TRUE) +
  theme_bw() +
  #geom_text(aes(x = rownames(dat)))+
  # scale_fill_manual(guide = "legend", labels=c("Triuret hydrolase", "Biuret hydrolase"), values = pal) +
  # scale_color_manual(guide="legend",labels=c("Triuret hydrolase", "Biuret hydrolase"),
  # values=pal) +
  theme(axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=20,face="bold",hjust=0)
  )
dev.off()

# PCA using FactoMineR package
# pca_2 <- PCA(findat, graph = FALSE)
# pca_2$var$coord
# 
# barplot(pca_2$eig[,1])


