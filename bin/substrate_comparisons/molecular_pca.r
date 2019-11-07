# Install packages
pacman::p_load("tidyverse", "ChemmineR", "readxl", "webchem", "Rcpi", 
               "recipes", "ggthemes", "caret", "earth")
# BiocManager::install("Rcpi", dependencies = c("Imports", "Enhances"))

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the raw data
rawdat <- read_csv("data/selected_molecular_properties_17pNPs.csv")

# Correlate with activity 
activity <- read_csv("output/substrate_comparisons/enzyme_activity_per_substrate.csv")
activity$abbrev

# Combine the data with activity
mergdat <- rawdat %>% 
  na.omit() %>%
  left_join(., activity, by = c("cmpnd_abbrev" = "abbrev"))
rownames(mergdat) <- mergdat$cmpnd_abbrev

# Trim irrelevant features
trimdat <- mergdat %>%
  dplyr::select(-status, -cmpnd_abbrev, -IUPAC, -SMILES) # %>%
  # as.data.frame()
rownames(trimdat)


# Principal components analysis in R
pca_dat <- prcomp(trimdat[,1:(ncol(trimdat)-1)], center = TRUE, scale. = TRUE)

# Split into training and testing
training_dat <- trimdat[!is.na(trimdat$activity),]
training_dat

testing_dat <- trimdat[is.na(trimdat$activity),]

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

# Put components axes on the same range

pca_rng <- extendrange(c(pca_test$PC1, pca_test$PC2))
ggplot(pca_test, aes(x = PC1, y = PC2)) +
  geom_point(alpha = .2, cex = 1.5) + 
  theme(legend.position = "top") +
  theme_bw() +
  xlim(pca_rng) + ylim(pca_rng) + 
  xlab("Principal Component 1") + ylab("Principal Component 2")


# Training a MARS model
mars_grid <- expand.grid(degree = 1:2, nprune = seq(2, 26, by = 2))



basic_rec <- recipe(activity ~ ., data = training_dat) %>%
  step_zv(all_predictors()) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  # step_pca(all_predictors()) %>%
  step_dummy(all_nominal())
  
ctrl <- trainControl(
  method = "cv", 
  # Save the assessment predictions from the best model
  savePredictions = "final",
  # Log the progress of the tuning process
  verboseIter = TRUE
)
training_dat

mars_mod <- train(
  data = training_dat,
  basic_rec,
  method = "earth",
  tuneGrid = mars_grid# ,
  # trControl = ctrl
)

mars_imp <- varImp(mars_mod)
ggplot(mars_imp, top = 20) + xlab("")
warnings()



rpart_mod <- train(
  data = training_dat,
  basic_rec,
  method = "rpart"#,
  #tuneLength = 10,
  #trControl = ctrl
)

par(xpd = NA) # Avoid clipping the text in some device
plot(rpart_mod$finalModel)
text(rpart_mod, digits = 3)

# Lasso mod,
# rpart_mod <- train(
#   data = training_dat,
#   basic_rec,
#   method = "enet",
#   alpha = 1,
#   trControl = trainControl(method = "cv", number = 10))
  
  #tuneLength = 10,
  #trControl = ctrl


knn_mod <- train(
  data = training_dat,
  basic_rec,
  method = "knn",
  tuneGrid = expand.grid(k = seq(1, 5, by = 1))) #,
  #trControl = trainControl(method = "cv", number = 5)
get_best_result(sim_knn_mod)


glm_mod <- train(
  data = training_dat,
  basic_rec,
  method = "glm",
  trControl = trainControl(method = "cv", number = 10))

# glm_mod <- train(
#   y = training_dat$activity,
#   x = training_dat[,2:(ncol(training_dat)-1)],
#   preProcess = c("center", "scale"),
#   method = "glmnet",
#   trControl = trainControl(method = "cv", number = 10))

glm_mod$recipe
glm_mod$results
important <- names(glm_mod$finalModel$coefficients)[!is.na(glm_mod$finalModel$coefficients)]
which_coeffs <- important[2:length(important)]
which_coeffs


library("PerformanceAnalytics")


coeff_dat <- training_dat %>%
  dplyr::select(which_coeffs)

grep("TPS", colnames(training_dat))
colnames(training_dat)
unscaled_dat <- training_dat %>%
  dplyr::select(c("activity", which_coeffs))

pdf("output/substrate_comparisons/molecular_descriptor_analytics_unscaled.pdf")
chart.Correlation(unscaled_dat, pch = 19, histogram = TRUE)
dev.off()

# find mean and sd column-wise of training data
trainMean <- apply(coeff_dat, 2, mean)
trainSd <- apply(coeff_dat,2,sd)


## centered AND scaled
norm2.testData <- sweep(sweep(coeff_dat, 2L, trainMean), 2, trainSd, "/")
norm2.testData

my_data <- cbind(norm2.testData, training_dat$activity)
colnames(my_data)[ncol(my_data)] <- 'activity'

pdf("output/substrate_comparisons/molecular_descriptor_analytics_scaled.pdf")
chart.Correlation(my_data, pch = 19, histogram = TRUE)
dev.off()

# GLmnet in caret 

# glm_mod <- train(
#   data = training_dat,
#   basic_rec,
#   method = "glmnet",
#   alpha = 1)
  # trControl = trainControl(method = "cv", number = 10))

,