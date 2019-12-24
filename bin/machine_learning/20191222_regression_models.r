# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
                "randomForest","gbm","nnet","e1071","svmpath","lars",
                "glmnet","svmpath")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the raw data
rawdat <- read_csv("data/machine_learning/20191222_1095_training_examples_pass1.csv")

# Only keep variables with nonzero variance
nozdat <- nearZeroVar(rawdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE]
which_rem

### INCLUDING SEQUENCE AND CHEMICAL FEATURES
dat <- rawdat %>%
  dplyr::mutate(id = paste0(org, "_", substrate)) %>%
  dplyr::select(-which_rem, -org, -substrate) %>%
  dplyr::select(id, contains("_"),  contains("PC"), activity)
dim(dat) #338 variables remaining with non-zero variance

# Set random seed 
set.seed(20192212)

# Split into test and training data
dat_split <- initial_split(dat, strata = "activity")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("id", "activity")]
x_test <- dat_test[,!colnames(dat_test) %in% c("id", "activity")]
y_train <- dat_train$activity
y_test <- dat_test$activity

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)

## First test a generalized linear model! # No penalization
glm_mod <- train(
  x = df_train,
  y = y_train,
  method = "glm",
  trControl = trainControl(method = "repeatedcv", 
                           repeats = 3, 
                           number = 10))

glm_mod$results

important <- names(glm_mod$finalModel$coefficients)[!is.na(glm_mod$finalModel$coefficients)]
which_coeffs <- important[2:length(important)]
glm_mod$finalModel$coefficients

# Now test elastic net (glmnet)
enet_mod <- train(
  x = df_train,
  y = y_train,
  method = "glmnet",
  trControl = trainControl(method = "repeatedcv", 
                           repeats = 3, 
                           number = 10))

enet_mod$results
important <- names(enet_mod$finalModel$coefficients)[!is.na(enet_mod$finalModel$coefficients)]
important
which_coeffs <- important[2:length(important)]

## Next test MARS
mars_grid <- expand.grid(degree = 1:2, nprune = seq(2, 26, by = 2))

mars_mod <- train(
  x = df_train,
  y = y_train,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = trainControl(method = "repeatedcv", 
                           repeats = 3, 
                           number = 10))
)

mars_imp <- varImp(mars_mod)
ggplot(mars_imp, top = 20) + xlab("")
mars_mod$results

