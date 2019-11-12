# Install packages
pacman::p_load("tidyverse", "ChemmineR", "readxl", "webchem", "Hmisc",
               "Rcpi", "recipes", "ggthemes", "caret", "earth")
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

# Correlation test
method <- "pearson"
x <- as.matrix(pca_dat$x, rownames.force = T)
x
y <- as.matrix(training_dat$activity, rownames.force = T)
colnames(y) <- "activity"
y
# x <- scale(x, center = T) # DO I NEED TO SCALE AND CENTER?

corrs <- list()
pvals <- list()
for(i in 1:ncol(x)){
    a <- x[,i,drop=F]
    b <- y
   # ll[[i]]$var1 <- colnames(x)[i]
    corrs[[i]] <- cor(a, b, use = "everything", method = method)
    pvals[[i]] <- cor.test(a, b, method = method)$p.value
}

findat <- data.frame(cbind(colnames(x), unlist(corrs), unlist(pvals)), stringsAsFactors = F)
colnames(findat) <- c("var", "corrs", "pvals")
findat$adj.pval.BH <- p.adjust(findat$pvals, method="BH") # Benjamin-Hochberg correction
# findat$adj.pval.boniferonni <- p.adjust(findat$pvals, method="bonferroni") #Bonferroni

which_sig <- findat[findat$pvals < 0.1,]
which_sig


### TRY IT WHERE THE COLUMNS ARE ADJUSTED
# x <- sweep(x, 2, colSums(x), '/')
x <- as.matrix(pca_dat$x, rownames.force = T)
# y <- y/colSums(y)
ylog <- log10(y)

colnames(y) <- "activity"
# x <- scale(x, center = T) # DO I NEED TO SCALE AND CENTER?

corrs <- list()
pvals <- list()
for(i in 1:ncol(x)){
  a <- x[,i,drop=F]
  b <- ylog[,1]
  # ll[[i]]$var1 <- colnames(x)[i]
  corrs[[i]] <- cor(a, b, use = "everything", method = method)
  pvals[[i]] <- cor.test(a, b, method = method)$p.value
}

findat.swept <- data.frame(cbind(colnames(x), unlist(corrs), unlist(pvals)), stringsAsFactors = F)
colnames(findat.swept) <- c("var", "corrs", "pvals")
findat.swept$adj.pval.BH <- p.adjust(findat.swept$pvals, method="BH")
as.numeric(findat.swept$pvals)# Benjamin-Hochberg correction
which_sig <- findat.swept[as.numeric(findat.swept$pvals) < 0.2,]
which_sig
which_corr <- findat.swept[as.numeric(findat.swept$corrs) > 0.25,]
which_corr

corrdat <- training_dat[,colnames(training_dat) %in% c("activity", which_corr$var)]
sigdat <- training_dat[,colnames(training_dat) %in% c("activity", which_sig$var)]


library("PerformanceAnalytics")
pdf("output/substrate_comparisons/molecular_descriptor_analytics_PC_all.pdf", width = 10, height = 10)
ch <- chart.Correlation(data.frame(cbind(x,y), stringsAsFactors = F), pch = 19, histogram = TRUE)
ch
dev.off()


