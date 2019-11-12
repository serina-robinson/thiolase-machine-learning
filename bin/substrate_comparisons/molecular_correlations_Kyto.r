# Install packages
pacman::p_load("tidyverse", "ChemmineR", "readxl", "webchem", "Hmisc",
               "Rcpi", "recipes", "ggthemes", "caret", "earth")
# BiocManager::install("Rcpi", dependencies = c("Imports", "Enhances"))

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the raw data
rawdat <- read_csv("data/selected_molecular_properties_17pNPs.csv")

# Correlate with all activity 
activity_all <- read_csv("output/substrate_comparisons/all_slopes_long_format.csv") 
activity_kyto <- activity_all %>%
  dplyr::filter(grepl("Kytococcus", org)) %>%
  dplyr::select(max_slope, cmpnd) %>%
  dplyr::rename(activity = max_slope)

# Combine the data with activity
mergdat <- rawdat %>% 
  na.omit() %>%
  left_join(., activity_kyto, by = c("cmpnd_abbrev" = "cmpnd"))
mergdat
rownames(mergdat) <- mergdat$cmpnd_abbrev

# Trim irrelevant features
trimdat <- mergdat %>%
  dplyr::select(-status, -cmpnd_abbrev, -IUPAC, -SMILES) # %>%
# as.data.frame()
trimdat[,ncol(trimdat)]


# Principal components analysis in R
# pca_dat <- prcomp(trimdat[,1:(ncol(trimdat)-1)], center = TRUE, scale. = TRUE)

# Split into training and testing
training_dat <- trimdat[!is.na(trimdat$activity),]
dim(training_dat)

testing_dat <- trimdat[is.na(trimdat$activity),]

# Correlation test
method <- "pearson"
feats <- training_dat[,1:(ncol(training_dat) - 1)]
nov <- feats[ - as.numeric(which(apply(feats, 2, var) == 0))] # remove columns with zero variance
x <- as.matrix(nov, rownames.force = T)
x
y <- as.matrix(training_dat$activity, rownames.force = T)
colnames(y) <- "activity"
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
x <- sweep(x, 2, colSums(x), '/')
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
which_sig <- findat.swept[as.numeric(findat.swept$pvals) < 0.1,]
which_sig
which_corr <- findat.swept[as.numeric(findat.swept$corrs) > 0.25,]
which_corr
corrdat <- training_dat[,colnames(training_dat) %in% c("activity", which_corr$var)]
corrdat
which_int <- c("MW", "ALogP", "AMR", "MolP", "nRotB", "TopoPSA", "nAtomLAC", "VABC")
intdat <- training_dat[,colnames(training_dat) %in% c("activity", which_int)]
intdat

###  Plot stuff
library("PerformanceAnalytics")

pdf("output/substrate_comparisons/molecular_descriptor_analytics_interesting.pdf", width = 10, height = 10)
ch <- chart.Correlation(intdat, pch = 19, histogram = TRUE)
ch
dev.off()

pdf("output/substrate_comparisons/molecular_descriptor_analytics_correlations.pdf", width = 10, height = 10)
ch <- chart.Correlation(corrdat, pch = 19, histogram = TRUE)
ch
dev.off()

pdf("output/substrate_comparisons/molecular_descriptor_analytics_significant.pdf", width = 100, height = 100)
ch <- chart.Correlation(corrdat, pch = 19, histogram = TRUE)
ch
dev.off()


### 

explore_pca <- recipe(activity ~ ., data = training_dat) %>%
  # step_BoxCox(all_predictors()) %>%
  # step_naomit(all_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  step_pca(all_predictors()) %>%
  step_dummy(all_nominal()) %>%
  prep(training = testing_dat, verbose = TRUE)

pca_test <- bake(explore_pca, new_data = testing_dat)
