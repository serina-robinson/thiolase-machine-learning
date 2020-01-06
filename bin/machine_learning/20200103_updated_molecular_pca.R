# Install packages
pacman::p_load("tidyverse", "ChemmineR", "readxl", "webchem", "Rcpi", 
               "recipes", "ggthemes", "caret", "earth", "factoextra", 
               "FactoMineR")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the raw data
rawdat <- read_csv("data/substrate_comparisons/15pNPs_159_selected_molecular_properties.csv") 

trimdat <- rawdat %>%
  dplyr::select(-cmpnd_abbrev, -IUPAC, -SMILES) # %>%
rownames(trimdat) <- rawdat$cmpnd_abbrev

# Only keep variables with nonzero variance
nozdat <- nearZeroVar(trimdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE]
which_rem

findat <- trimdat %>%
  dplyr::select(-which_rem)
head(findat)

# Principal components analysis in R
set.seed(20200104)
pca_dat <- prcomp(findat, center = TRUE, scale. = TRUE)
pca_abs <- abs(pca_dat$rotation) #Pull the absolute values of the rotations
head(pca_dat$rotation)

screeplot(pca_dat, npcs = 15, type = "lines")

pdf("output/substrate_comparisons/molecular_descriptor_PC_screeplot_variances.pdf")
screeplot(pca_dat, npcs = 15, type = "lines")
dev.off()

imp <- summary(pca_dat)[6][[1]] 
sum(imp[2,1:7]) # First 7 PCs explain over 92% of the total variances

## Try grab the first 7 PCs
# pca_abs[,4][order(pca_abs[,4])]
# pca_abs[order(pca_abs[,3], decreasing = T)]

res <- lapply(1:7, function(x) { rownames(pca_abs)[order(pca_abs[,x], decreasing = T)] })
res_load <- lapply(1:7, function(x) { pca_abs[,x][order(pca_abs[,x], decreasing = T)] })
res_load

names(res) <- paste0("PC", 1:7)
names(res_load) <- paste0("PC_loading", 1:7)

res_df <- bind_rows(res)
res_load_df <- bind_rows(res_load)

findf <- bind_cols(res_df, res_load_df) %>%
  select(ends_with("1"), ends_with("2"), ends_with("3"), ends_with("4"), ends_with("5"), ends_with("6"), ends_with("7"))
findf
write_csv(findf, "data/machine_learning/PC7_loadings.csv")


write_csv(data.frame(cbind(rownames(pca_dat$x), pca_dat$x[,1:7]), stringsAsFactors = F), 
          "data/machine_learning/PC7_molecular_descriptors.csv")


dat <- data.frame(pca_dat$x[,1:2], stringsAsFactors = F)
dat


pdf(paste0("output/substrate_comparisons/molecular_descriptor_pca_dat.pdf"), width = 5, height = 5)
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


