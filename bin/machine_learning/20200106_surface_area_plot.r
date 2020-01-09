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
  dplyr::select(substrate, nB, MW,  AROMATIC, O, N,  MLogP, nRotB, Cl, VABC, nAtomLAC, TopoPSA, HybRatio)
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
  dplyr::select(org, substrate, sub_num, avg_sub, volume_sa, area_sa, contains("PC"), HybRatio, VABC, nAtomLAC, MLogP, TopoPSA, O, MW, nRotB, avg_activity, activity, n) %>%
  # dplyr::select(org, volume_sa, area_sa, avg_activity, n) %>%
  dplyr::rename(`Number of pNP substrates accepted` = n) %>%
  dplyr::rename(`Average activity` = avg_activity) %>%
  group_by(substrate) %>%
  # dplyr::mutate(num_sub = paste0(rep(1:14, each = n), substrate)) %>%
  distinct() 

tma <- combdat %>%
  dplyr::filter(substrate == "TMA")


dimethyl <- combdat %>%
  dplyr::filter(substrate == "dimethyl")

my.formula <- y ~ x
pdf("output/machine_learning/solvent_pocket_volume_surface_area_plot_TMA.pdf", width = 8, height = 8)
ggplot(data = tma, aes(x = area_sa, y = volume_sa)) +
  geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  geom_point(aes(color = activity,
                 size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  ggrepel::geom_text_repel(data = subset(tma, activity > 1.25), segment.size = 0.2,
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


pl1 <- ggplot(data = combdat, aes(x = activity, y = sub_num, group = sub_num)) +
  #geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  #geom_point(aes(color = VABC, alpha = 0.7)) +
  #size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  scale_fill_viridis(option = "plasma") +
  geom_density_ridges_gradient(aes(color = HybRatio, fill = HybRatio),
                               jittered_points = TRUE,
                               quantile_lines = TRUE, scale = 0.6, alpha = 0.7,
                               vline_size = 1, vline_color = "gray40",
                               point_size = 0.4, point_alpha = 1,
                               position = position_raincloud(adjust_vlines = TRUE)) +
  xlab("Enzyme activity") +
  ylab("Substrate")
pl1

pdf("output/machine_learning/solvent_pocket_volume_surface_area_plot_HybRatio.pdf", width = 8, height = 8)
ggplot(data = combdat, aes(y = area_sa, x = HybRatio)) +
  geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  geom_point(aes(color = activity,
                 size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  # ggrepel::geom_text_repel(data = subset(combdat, activity > 1.25), segment.size = 0.2,
  #                          nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
  #                          segment.color = "grey50", #direction = "x",
  #                          aes(label = org), fontface = "italic") +
  ylab("Solvent-accessible surface area") +
  xlab("Hybridization ratio") +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) 
#scale_color_manual(values = RColorBrewer::brewer.pal(8, "Spectral"))
dev.off()

pdf("output/machine_learning/solvent_pocket_volume_surface_area_plot_dimethyl.pdf", width = 8, height = 8)
ggplot(data = dimethyl, aes(x = area_sa, y = volume_sa)) +
  geom_smooth(method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  geom_point(aes(color = activity,
                 size = activity), alpha = 0.7) +
  theme_pubr() +
  scale_color_viridis(option = "plasma") +
  ggrepel::geom_text_repel(data = subset(dimethyl, activity > 1.25), segment.size = 0.2,
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