# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest","gbm","nnet","e1071","lars",
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
  dplyr::mutate(is_active = as.factor(case_when(activity == 0 ~ "N",
                                      TRUE ~ "Y"))) %>%
  dplyr::select(-which_rem, -org, -substrate) %>%
  dplyr::select(id, contains("_"),  contains("PC"), is_active)
dim(dat) #338 variables remaining with non-zero variance
dat$is_active

# Set random seed 
set.seed(20192212)

# Split into test and training data
dat_split <- initial_split(dat, strata = "is_active")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("id", "is_active")]
x_test <- dat_test[,!colnames(dat_test) %in% c("id", "is_active")]
y_train <- dat_train$is_active
y_test <- dat_test$is_active

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$id)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$id)

# Quick  test

rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                          repeats = 3,
                          verboseIter = T, classProbs = T,
                          savePredictions = "final"),
  # tuneGrid = tgrid,
  num.trees = 500,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  importance = "permutation") 

# Confusion matrix
getTrainPerf(rf) # Training set accuracy is 80% for random forest

approx_roc_curve <- function(x, label) {
  x %>%
    pluck("pred") %>%
    roc_curve(obs, y_train) %>%
    mutate(model = label)
}

approx_roc_curve(rf, "Random Forest") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)


# Try prediction
rf_ml <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = as.character(rf$bestTune$splitrule),
                mtry = rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
                importance = "permutation")

rf_pred <- predict(rf, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$is_active))
cm_rf

sink("output/machine_learning/random_forest_binary_classification_results.txt")
cm_rf
sink()

vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T),
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:25)
vimp

pdf("output/machine_learning/random_forest_binary_classification_vimp.pdf", width = 6, height = 6)
ggplot(data = vimp,
       aes(x=reorder(variable,importance), y=importance, fill=importance))+
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
dev.off()

###Naive Bayes method
nb_grid <- expand.grid(usekernel = TRUE, fL = 0, adjust = 1)

nb <- train(
  x = df_train, 
  y = y_train,
  method = "nb",
  tuneGrid = nb_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, 
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"))

# Confusion matrix
getTrainPerf(nb) # Only 59% training accuracy with NB
nb$results


nb$results
nb_pred <- predict(nb, newdata = form_test)
warnings()
cm_nb <- confusionMatrix(nb_pred, as.factor(dat_test$is_active))

sink("output/machine_learning/Naive_Bayes_binary_classification_results.txt")
cm_nb
sink()

## Make a heatmap of confusion matrix results
cm_list <- list(cm_nb$table)
names(cm_list) <- c("nb_aa")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("output/machine_learning/nb_binary_classification_matrix_20191224.jpeg", height=7, width=7, units='in')

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
                           savePredictions = "final"),
  preProcess = c("center", "scale")
)

# Confusion matrix
getTrainPerf(nnet_mod)


# size decay
#  10 0.5

# approx_roc_curve(nnet_ml, "Neural network") %>%
#   ggplot(aes(x = 1 - specificity, y = sensitivity)) +
#   geom_path()  +
#   geom_abline(col = "red", alpha = .5)

nnet_pred <- predict(nnet_mod, newdata = form_test)
cm_nnet <- confusionMatrix(nnet_pred, as.factor(dat_test$is_active))
cm_nnet

sink("output/machine_learning/Neural_network_binary_classification_results.txt")
cm_nnet
sink()
# dtl_feat_select <- data.frame(round(cm_nnet$byClass[,colnames(cm_nnet$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))

## Make a heatmap of confusion matrix results
cm_list <- list(cm_nnet$table)
names(cm_list) <- c("nnet_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("output/machine_learning/nnet_binary_classification_matrix_20191224.jpeg", height=7, width=7, units='in')

###==================================================================###
###SVM Radial basis kernel
# Fitting sigma = 0.000646, C = 1 on full training set
grid <- expand.grid(.C = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1),
                    .sigma = c(0.1, 0.25,0.5,0.75,1))

svm_radial <- train(
  x = df_train,
  y = y_train,
  method = "svmRadial",
  tuneGrid = grid,
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  preProcess = c("center", "scale"),
  verbose = TRUE)

# Confusion matrix
getTrainPerf(svm_radial)
svm_radial$bestTune

approx_roc_curve(svm_radial, "SVM Radial") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

pdf("output/machine_learning/svm_radial_feat_select_var_imp_binary_classification.pdf", width = 6, height = 6)
ggplot(data = vimp, 
       aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
dev.off()

cm_svm_radial_pred <- predict(svm_radial, newdata = form_test) # FAILS FOR SOME REASON
cm_svm_radial <- confusionMatrix(svm_radial_pred, as.factor(dat_test$clf))

dtl_feat_select <- data.frame(round(cm_svm_radial$byClass[,colnames(cm_svm_radial$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))

## Make a heatmap of confusion matrix results
cm_list <- list(cm_svm_radial$table)
names(cm_list) <- c("svm_radial_aa")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("output/machine_learning/svm_radial_conf_matrices_binary_classification_20191228.jpeg", height=7, width=7, units='in')

###==================================================================###
###SVM Linear kernel
# Fitting sigma = 0.000646, C = 1 on full training set


# grid <- expand.grid(C = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1),
#                     gamma=c(0.1, 0.25,0.5,0.75,1))
# svm_linear <- train(
#   x = df_train,
#   y = y_train,
#   tuneGrid= grid,
#   tuneLength = 10,
#   method = "svmLinear",
#   trControl = trainControl(method = "repeatedcv", 
#                            number = 10, repeats = 3,
#                            verboseIter = T, classProbs = T,
#                            savePredictions = "final"),
#   preProcess = c("center", "scale"),
#   verbose = TRUE)
# 
# # Confusion matrix
# getTrainPerf(svm_linear)
# 
# approx_roc_curve(svm_linear, "SVM Linear") %>%
#   ggplot(aes(x = 1 - specificity, y = sensitivity)) +
#   geom_path()  +
#   geom_abline(col = "red", alpha = .5)
# 
# pdf("output/machine_learning/svm_linear_feat_select_var_imp_binary_classification_20191228.pdf", width = 6, height = 6)
# ggplot(data = vimp, 
#        aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
#   geom_bar(stat="identity", position="dodge")+ coord_flip()+
#   ylab("Variable Importance")+
#   xlab("")+
#   guides(fill=F)+
#   scale_fill_gradient(low="red", high="blue")
# dev.off()
# 
# cm_svm_linear_pred <- predict(svm_linear, newdata = form_test)
# cm_svm_linear <- confusionMatrix(svm_linear_pred, as.factor(dat_test$is_active))
# 
# 
# dtl_feat_select <- data.frame(round(cm_svm_linear$byClass[,colnames(cm_svm_linear$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
# 
# ## Make a heatmap of confusion matrix results
# cm_list <- list(cm_svm_linear$table)
# names(cm_list) <- c("svm_linear_aa")
# 
# pllist <- list()
# for(i in 1:length(cm_list)) {
#   pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
#     geom_tile(aes(fill = Freq)) + 
#     geom_text(aes(label = round(Freq, 1))) +
#     scale_fill_gradient(low = "white", high = "red") +
#     ggtitle(names(cm_list)[i]) +
#     theme(axis.text.x = element_text(angle = 90), 
#           axis.title.x = element_blank())
# }
# 
# pllist[[1]]
# ggsave("output/machine_learning/svm_linear_binary_classification_matrix_20191228.jpeg", height=7, width=7, units='in')
# 
# ##  Compare all confusion matrices
# cm_list <- list(cm_rf$table,
#                 cm_nb$table,
#                 cm_nnet$table,
#                 cm_svm_radial$table,
#                 cm_svm_linear$table)
# names(cm_list) <- c("rf", "nb", "nnet", "svm_radial", "svm_linear")
# 
# pllist <- list()
# for(i in 1:length(cm_list)) {
#   pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
#     geom_tile(aes(fill = Freq)) + 
#     geom_text(aes(label = round(Freq, 1))) +
#     scale_fill_gradient(low = "white", high = "red") +
#     ggtitle(names(cm_list)[i]) +
#     theme(axis.text.x = element_text(angle = 90), 
#           axis.title.x = element_blank())
#   #axis.title.y = element_blank())
# }
# 
# pllist[[1]]
# 
# plot_grid(pllist[[1]],
#           pllist[[2]],
#           pllist[[3]],
#           pllist[[4]],
#           pllist[[5]])
# ggsave("output/machine_learning/model_confusion_matrices.jpeg", height=6, width=10, units='in')

