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

### INCLUDING BOTH SEQUENCE AND CHEMICAL FEATURES
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
                           number = 10, 
                           savePredictions = ))

glm_mod$results # RMSE 0.679 is pretty bad

important <- names(glm_mod$finalModel$coefficients)[!is.na(glm_mod$finalModel$coefficients)]
which_coeffs <- important[2:length(important)]
glm_mod$finalModel$coefficients
glm_mod$pred

# Training and test performance
getTrainPerf(glm_mod)
glm_pred <- predict(glm_mod, newdata = form_test)


# Now test elastic net (glmnet)
enet_mod <- train(
  x = df_train,
  y = y_train,
  method = "glmnet",
  trControl = trainControl(method = "repeatedcv", 
                           repeats = 3, 
                           number = 10))

enet_mod$results # best RMSE is 0.7
important <- names(enet_mod$finalModel$coefficients)[!is.na(enet_mod$finalModel$coefficients)]
important
which_coeffs <- important[2:length(important)]

# Confusion matrix
getTrainPerf(enet_mod)
enet_pred <- predict(enet_mod, newdata = form_test)



# Random forest regression
rf_mod <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           savePredictions = "final"),
  # tuneGrid = tgrid,
  num.trees = 1000,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  importance = "permutation") # BEST R MSE is 55.7%

# Confusion matrix
getTrainPerf(rf_mod)
enet_pred <- predict(enet_mod, newdata = form_test)

# NOT RUN
# approx_roc_curve <- function(x, label) {
#   x %>%
#     pluck("pred") %>%
#     roc_curve(obs, y_train) %>%
#     mutate(model = label)
# }
# 
# approx_roc_curve(rf, "Random Forest") %>%
#   ggplot(aes(x = 1 - specificity, y = sensitivity)) +
#   geom_path()  +
#   geom_abline(col = "red", alpha = .5)
# 
# # Confusion matrix
# getTrainPerf(rf) # RMSE is 55.7%
# 
# # Try prediction
# rf_ml <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = as.character(rf$bestTune$splitrule), 
#                 mtry = rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
#                 importance = "permutation")
# 
# rf_pred <- predict(rf, newdata = form_test)
# cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$clf))
# 
# sink("output/machine_learning/random_forest_regression_features.txt")
# cm_rf
# sink()
# 
# vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T), 
#                          names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
#   dplyr::rename(importance = X1,
#                 variable = X2) %>%
#   mutate(importance = as.numeric(importance)) %>%
#   dplyr::slice(1:25)
# vimp
# 
# pdf("output/machine_learning/random_forest_regression.pdf", width = 6, height = 6)
# ggplot(data = vimp, 
#        aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
#   geom_bar(stat="identity", position="dodge")+ coord_flip()+
#   ylab("Variable Importance")+
#   xlab("")+
#   guides(fill=F)+
#   scale_fill_gradient(low="red", high="blue")
# dev.off()
# 
# 
# 
# obs_pred_plot <- function(x, dat, cutoff = 25, ...) {
#   pred_dat <- x %>%
#     add_columns(dat, model, year) %>%
#     mutate(residuals = obs - pred) 
#   ggplot(pred_dat, aes(x = pred, y = obs)) +
#     geom_abline(col = "green", alpha = .5) + 
#     geom_point(alpha = .3) + 
#     geom_smooth(
#       se = FALSE, col = "red", 
#       lty = 2, lwd = .25, alpha = .5
#     ) + 
#     geom_text_repel(
#       data = dplyr::filter(pred_dat, abs(residuals) > cutoff),
#       aes(label = plot_label),
#       segment.color = "grey50"
#     )
# }

## Next test MARS
mars_grid <- expand.grid(degree = 1:2, nprune = seq(2, 26, by = 2))

mars_mod <- train(
  x = df_train,
  y = y_train,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = trainControl(method = "repeatedcv", 
                           repeats = 3, 
                           number = 10,
                           savePredictions = "final"))
)

mars_imp <- varImp(mars_mod)
ggplot(mars_imp, top = 20) + xlab("")
mars_mod$results
ggplot(mars_mod$pred, aes(x = obs, y = pred)) +
  geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", 
              lty = 2, lwd = 1, alpha = .5)

obs_pred_plot <- function(x, dat, cutoff = 25, ...) {
  pred_dat <- x %>%
    add_columns(dat, model, year) %>%
    mutate(residuals = obs - pred) 
  ggplot(pred_dat, aes(x = pred, y = obs)) +
    geom_abline(col = "green", alpha = .5) + 
    geom_point(alpha = .3) + 
    geom_smooth(
      se = FALSE, col = "red", 
      lty = 2, lwd = .25, alpha = .5
    ) + 
    geom_text_repel(
      data = dplyr::filter(pred_dat, abs(residuals) > cutoff),
      aes(label = plot_label),
      segment.color = "grey50"
    )
}

###Neural networks using nnet
# Currently throws an error
nnet_grid <- expand.grid(.decay = c(0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7),
                         .size = c(3, 5, 10, 20))

nnet_mod <- train(
  x = df_train,
  y = y_train,
  method = "nnet",
  MaxNWts = 1000000000,
  # size = 1,
  tuneGrid = nnet_grid,
  trControl = trainControl(method = "repeatedcv", repeats = 3,
                           number = 10,
                           verboseIter = T,
                           savePredictions = "final"))

# Confusion matrix
getTrainPerf(nnet_mod)

# size decay
#  10 0.5

# approx_roc_curve(nnet_ml, "Neural network") %>%
#   ggplot(aes(x = 1 - specificity, y = sensitivity)) +
#   geom_path()  +
#   geom_abline(col = "red", alpha = .5)

nnet_pred <- predict(nnet_mod, newdata = form_test)
cm_nnet <- confusionMatrix(nnet_pred, as.factor(dat_test$activity))
cm_nnet

sink("output/machine_learning/Neural_network_regression_results.txt")
cm_nnet
sink()


