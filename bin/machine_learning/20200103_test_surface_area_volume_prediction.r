# Read in the raw data
# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest","gbm","nnet","e1071","svmpath","lars",
               "glmnet","svmpath", "readxl", "ggridges", "viridis", "ggpubr", "ggpmisc")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the molecular features
molec_fts <- read_csv("data/machine_learning/PC7_molecular_descriptors.csv") %>%
  dplyr::rename(substrate = V1)
molec_fts$substrate <- gsub(" ", "\\.", molec_fts$substrate)
molec_fts$substrate[molec_fts$substrate == "oxidazole"] <- "oxadiazole"

loadings <- read_csv("output/substrate_comparisons/20191218_relative_PC_loadings.csv")
chem_descr <- read_csv("data/substrate_comparisons/15pNPs_159_selected_molecular_properties.csv") %>%
  dplyr::mutate(substrate = gsub(" ", "\\.", cmpnd_abbrev)) %>%
  dplyr::select(substrate, nB, MW,  AROMATIC, O, N,  MLogP, nRotB, Cl, VABC, nAtomLAC, TopoPSA)
chem_descr$substrate[chem_descr$substrate == "oxidazole"] <- "oxadiazole"

# Read in the activity data
activity <- read_csv("data/machine_learning/20191218_all_cmpnds_avg_log_slopes_for_modeling.csv")

# Read in the cavity volumes
cavity_fts <- read_excel("data/caver_models/CastP_volume_rawdata.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::mutate(genus = gsub("\\.pdb", "", word(filename, sep = "_", 4))) %>%
  dplyr::mutate(species = gsub("\\.pdb", "", word(filename, sep = "_", 5))) %>%
  dplyr::mutate(org = paste0(genus, " ", species)) 

# Now merge everything...
comb <- activity %>%
  dplyr::left_join(., molec_fts, by = "substrate") %>%
  dplyr::left_join(., cavity_fts, by = "org") %>%
  dplyr::left_join(., chem_descr, by = "substrate") %>%
  dplyr::mutate(id = paste0(org, "_", substrate)) %>%
  dplyr::mutate(is_active = as.factor(case_when(activity == 0 ~ "N",
                                                TRUE ~ "Y"))) 

comb$area_sa <- as.numeric(comb$area_sa)
comb$volume_sa <- as.numeric(comb$volume_sa)

avgdat <- comb %>%
  dplyr::filter(activity != 0) %>%
  group_by(org) %>%
  summarise_at(vars(activity), mean) %>%
  dplyr::rename(avg_activity = activity)

avgsub <- comb %>%
  dplyr::filter(activity != 0) %>%
  group_by(substrate) %>%
  summarise_at(vars(activity), median) %>%
  dplyr::rename(avg_sub = activity) %>%
  arrange(desc(avg_sub)) %>%
  dplyr::mutate(sub_num = paste0(letters[1:14], "_", substrate))
avgsub

combdat <- comb %>%
  right_join(avgsub, ., by = "substrate") %>%
  left_join(., avgdat, by = "org") %>%
  dplyr::filter(activity != 0) %>%
  add_count(org) %>%
  #add_count(substrate) %>%
  dplyr::filter(to_keep == "Y") %>%
  dplyr::select(org, substrate, sub_num, avg_sub, volume_sa, area_sa, contains("PC"), VABC, nAtomLAC, MLogP, TopoPSA, O, MW, nRotB, avg_activity, activity, n) %>%
  # dplyr::select(org, volume_sa, area_sa, avg_activity, n) %>%
  dplyr::rename(`Number of pNP substrates accepted` = n) %>%
  dplyr::rename(`Average activity` = avg_activity) %>%
  group_by(substrate) %>%
  # dplyr::mutate(num_sub = paste0(rep(1:14, each = n), substrate)) %>%
  distinct() 



head(combdat)

# Now check for missing data and duplicates
dedup <- comb[complete.cases(comb),] 
dat <- dedup[!duplicated(comb),]

# For further analysis, remove errors in protein cavity calculations
outl <- dat %>% 
  dplyr::filter(to_keep == "Y")
hist(outl$area_sa)
hist(outl$volume_sa)

subdat <- outl[outl$is_active == "Y",]

# Plot the volume and SA vs. activity
# Molecular descriptors vs. activity

#combdat$substrate <- as.factor(combdat$substrate)

# levs <- as.character(unique(combdat$substrate))
# levels(combdat$substrate) <- rev(levs)

pl1 <- ggplot(data = combdat, aes(x = activity, y = sub_num, group = sub_num)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  scale_fill_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = VABC, fill = VABC),
                               jittered_points = TRUE,
                               quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
                               vline_size = 1, vline_color = "gray40",
                               point_size = 0.4, point_alpha = 1,
                               position = position_raincloud(adjust_vlines = TRUE)) +
  xlab("Enzyme activity") +
  ylab("Substrate")

pl1

pl2 <- ggplot(data = combdat, aes(x = activity, y = substrate, group = substrate)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  scale_fill_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = O, fill = O),
                               jittered_points = TRUE,
                               quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
                               vline_size = 1, vline_color = "gray40",
                               point_size = 0.4, point_alpha = 1,
                               position = position_raincloud(adjust_vlines = TRUE)) +
  xlab("Enzyme activity") +
  ylab("Substrate")

pl2

pl3 <- ggplot(data = combdat, aes(x = activity, y = substrate, group = substrate)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  scale_fill_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = MW, fill = MW),
                               jittered_points = TRUE,
                               quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
                               vline_size = 1, vline_color = "gray40",
                               point_size = 0.4, point_alpha = 1,
                               position = position_raincloud(adjust_vlines = TRUE)) +
  xlab("Enzyme activity") +
  ylab("Substrate")

pl3

pl5 <- ggplot(data = combdat, aes(x = activity, y = substrate, group = substrate)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  scale_fill_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = TopoPSA, fill = TopoPSA),
                               jittered_points = TRUE,
                               quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
                               vline_size = 1, vline_color = "gray40",
                               point_size = 0.4, point_alpha = 1,
                               position = position_raincloud(adjust_vlines = TRUE)) +
  xlab("Enzyme activity") +
  ylab("Substrate")

pl5

pl6 <- ggplot(data = combdat, aes(x = activity, y = substrate, group = substrate)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  scale_fill_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = MLogP, fill = MLogP),
                               jittered_points = TRUE,
                               quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
                               vline_size = 1, vline_color = "gray40",
                               point_size = 0.4, point_alpha = 1,
                               position = position_raincloud(adjust_vlines = TRUE)) +
  xlab("Enzyme activity") +
  ylab("Substrate")

pl6

pl7 <- ggplot(data = combdat, aes(x = activity, y = substrate, group = substrate)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  scale_fill_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = nAtomLAC, fill = nAtomLAC),
                               jittered_points = TRUE,
                               quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
                               vline_size = 1, vline_color = "gray40",
                               point_size = 0.4, point_alpha = 1,
                               position = position_raincloud(adjust_vlines = TRUE)) +
  xlab("Enzyme activity") +
  ylab("Substrate")

pl7

pl4 <- ggplot(data = combdat, aes(x = activity, y = substrate, group = substrate)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  scale_fill_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = nRotB, fill = nRotB),
                               jittered_points = TRUE,
                               quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
                               vline_size = 1, vline_color = "gray40",
                               point_size = 0.4, point_alpha = 1,
                               position = position_raincloud(adjust_vlines = TRUE)) +
  xlab("Enzyme activity") +
  ylab("Substrate")

pl4

pl8 <- ggplot(data = combdat, aes(x = activity, y = sub_num, group = sub_num)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  scale_fill_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = PC2, fill = PC2),
                               jittered_points = TRUE,
                               quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
                               vline_size = 1, vline_color = "gray40",
                               point_size = 0.4, point_alpha = 1,
                               position = position_raincloud(adjust_vlines = TRUE)) +
  xlab("Enzyme activity") +
  ylab("Substrate")

pl8

pl8 <- ggplot(data = combdat, aes(x = activity, y = sub_num, group = sub_num)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  scale_fill_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = PC2, fill = PC2),
                               jittered_points = TRUE,
                               quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
                               vline_size = 1, vline_color = "gray40",
                               point_size = 0.4, point_alpha = 1,
                               position = position_raincloud(adjust_vlines = TRUE)) +
  xlab("Enzyme activity") +
  ylab("Substrate")

pdf("output/PC2_ggridges.pdf")
pl8
dev.off()

pdf("output/chemical_descriptors_with_activity_distribution_ggridges.pdf", width = 25, height = 10)
plot_grid(pl1, pl3, pl4, pl5, pl6, pl7, ncol = 6)
dev.off()

pdf("output/machine_learning/solvent_pocket_volume_surface_area_geom_ridges_plot.pdf", width = 8, height = 15)
ggplot(data = combdat, aes(x = activity, y = substrate)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
                 #size = activity), alpha = 0.7) +
  theme_pubr() +
  #scale_color_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = VABC, fill = VABC),
    jittered_points = TRUE,
    quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
    vline_size = 1, vline_color = "red", vline_alpha = 1,
    point_size = 0.4, point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE)) 
dev.off()


pdf("output/machine_learning/solvent_pocket_volume_surface_area_facet_grid_plot.pdf", width = 8, height = 15)
ggplot(data = combdat, aes(x = area_sa, y = activity)) +
  geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  geom_density_ridges_gradient(
    jittered_points = TRUE,
    quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
    vline_size = 1, vline_color = "red", vline_alpha = 1,
    point_size = 0.4, point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE)) +
  ggrepel::geom_text_repel(data = subset(combdat, activity > 1.25), segment.size = 0.2,
                           nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
                           segment.color = "grey50", #direction = "x",
                           aes(label = org), fontface = "italic") +
  ylab("Solvent-accessible pocket volume") +
  xlab("Solvent-accessible surface area") +
  facet_grid(as.factor(combdat$substrate))
dev.off()

combdat <- comb %>%
  left_join(., avgdat, by = "org") %>%
  dplyr::filter(activity != 0) %>%
  add_count(org) %>%
  dplyr::filter(to_keep == "Y") %>%
  # dplyr::select(org, substrate, volume_sa, area_sa, VABC, avg_activity, activity, n) %>%
  dplyr::select(org, volume_sa, area_sa, avg_activity, n) %>%
  dplyr::rename(`Number of pNP substrates accepted` = n) %>%
  dplyr::rename(`Average activity` = avg_activity) %>%
  distinct()

my.formula <- y ~ x
pdf("output/machine_learning/solvent_pocket_volume_surface_area_plot.pdf", width = 8, height = 8)
ggplot(data = combdat, aes(x = area_sa, y = volume_sa)) +
  geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  geom_point(aes(color = `Number of pNP substrates accepted`,
                 size = `Average activity`), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  ggrepel::geom_text_repel(data = subset(combdat, `Average activity` > 1.25), segment.size = 0.2,
                           nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
                           segment.color = "grey50", #direction = "x",
                           aes(label = org), fontface = "italic") +
  ylab("Solvent-accessible pocket volume") +
  xlab("Solvent-accessible surface area") +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) 
#scale_color_manual(values = RColorBrewer::brewer.pal(8, "Spectral"))
dev.off()



pdf("output/machine_learning/solvent_pocket_volume_surface_area_plot.pdf", width = 8, height = 8)
ggplot(data = combdat, aes(x = area_sa, y = volume_sa)) +
  geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  geom_point(aes(color = `Number of pNP substrates accepted`,
             size = `Average activity`), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  ggrepel::geom_text_repel(data = subset(combdat, `Average activity` > 1.25), segment.size = 0.2,
                           nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
                           segment.color = "grey50", #direction = "x",
                           aes(label = org), fontface = "italic") +
  ylab("Solvent-accessible pocket volume") +
  xlab("Solvent-accessible surface area") +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) 
  #scale_color_manual(values = RColorBrewer::brewer.pal(8, "Spectral"))
dev.off()


ggplot(data = subdat,
       aes(x = volume_sa, y = activity)) +
  geom_point() +
  facet_grid(as.factor(subdat$substrate))

ggplot(data = subdat, aes(x = area_sa, y = activity)) +
  geom_point()

ggplot(data = subdat, aes(x = area_sa, y = volume_sa)) +
  geom_point()

ggplot(subdat, aes(x = volume_sa, y = substrate, fill = ..x..)) +
  geom_density_ridges_gradient(
    jittered_points = TRUE, 
    quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
    vline_size = 1, vline_color = "red", vline_alpha = 1,
    point_size = 0.4, point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE)) +
  scale_fill_viridis("Pocket volume") +
  theme_pubr() +
  xlab("Pocket solvent-accessible volume")

ggplot(subdat, aes(x = area_sa, y = substrate, fill = ..x..)) +
  #geom_density_ridges_gradient(scale = 0.8, rel_min_height = 0.01) +
  geom_density_ridges_gradient(
    jittered_points = TRUE, 
    quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
    vline_size = 1, vline_color = "red", vline_alpha = 1,
    point_size = 0.4, point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE)) +
  scale_fill_viridis("Pocket surface area") +
  theme_pubr() +
  xlab("Pocket solvent-accessible surface area")

ggplot(subdat, aes(x = activity, y = substrate, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 0.8, rel_min_height = 0.01) +
  scale_fill_viridis("Activity") +
  theme_pubr() +
  xlab("Activity")

ggplot(subdat, aes(x = activity, y = substrate, fill = ..x..)) +
  geom_density_ridges_gradient(
    jittered_points = TRUE, 
    quantile_lines = TRUE, scale = 0.9, alpha = 0.7,
    vline_size = 1, vline_color = "red", vline_alpha = 1,
    point_size = 0.4, point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE)) +
  scale_fill_viridis("Activity") +
  theme_pubr()

# Can compute the average surface area/volume of those that are active and compare to those that are inactive


ggplot(data = outl, aes(x = area_sa, y = volume_sa)) +
  geom_point(color = ifelse(outl$is_active == "Y", "maroon", "blue")) +
  facet_grid(as.factor(outl$substrate))

ggplot(data = outl, aes(x = area_sa, y = volume_sa)) +
  geom_point(color = ifelse(outl$is_active == "Y", "maroon", "blue"))
  #facet_grid(as.factor(outl$substrate))


# 


combdat <- outl %>%
  left_join(., avgdat, by = "org") %>%
  dplyr::filter(activity != 0) %>%
  add_count(org) %>%
  dplyr::select(org, volume_sa, area_sa, avg_activity, n) %>%
  distinct()

combdat$
ggplot(data = combdat, aes(x = area_sa, y = volume_sa)) +
  geom_point(alpha = 0.4, #color = as.factor(n),
             #color = ifelse(combdat$avg_activity > 1.5, "maroon", "blue"),
             size = combdat$avg_activity * 5) +
  theme_pubr() +
  ggrepel::geom_text_repel(data = subset(combdat, avg_activity > 1.5),
                           aes(label = org), color = "maroon")



ggplot(data = combdat, 
       aes(x = as.factor(n), y = area_sa)) +
  geom_point(aes(size = avg_activity),
             alpha = 0.4, 
             color = ifelse(combdat$avg_activity > 1.5, "maroon", "blue")) +
  theme_pubr() +
  ggrepel::geom_text_repel(data = subset(combdat, avg_activity > 1.5),
                           aes(label = org), color = "maroon")

ggplot(data = combdat, aes(x = area_sa, y = volume_sa)) +
  geom_point(alpha = 0.4, 
             color = ifelse(combdat$n >= 12, "maroon", "blue"),
             aes(size = avg_activity)) +
  theme_pubr() +
  ggrepel::geom_text_repel(data = subset(combdat, n >= 12),
                           aes(label = org), color = "maroon")






# Only keep variables with nonzero variance
nozdat <- nearZeroVar(dat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE]

# Set random seed 
set.seed(20200103)

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
table(dat$is_active)

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
  preProcess = c("center", "scale"),
  importance = "permutation") 

# Confusion matrix
getTrainPerf(rf) # Training set accuracy is 82% for random forest

# Try prediction
rf$bestTune$splitrule
rf$bestTune$mtry
ncol(form_train)

rf_ml <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini", #splitrule = as.character(rf$bestTune$splitrule),
                mtry = sqrt(ncol(form_train)), #rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
                importance = "permutation")

rf_pred <- predict(rf, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$is_active))
cm_rf # 80.22% test set accuracy

sink("output/machine_learning/random_forest_cavity_volume_binary_classification_results_20200103.txt")
cm_rf
sink()

vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T),
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:15)
vimp

ggplot(data = vimp,
       aes(x=reorder(variable,importance), y=importance, fill=importance))+
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")

