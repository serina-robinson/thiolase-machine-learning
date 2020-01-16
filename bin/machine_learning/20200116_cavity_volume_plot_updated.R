# Read in the raw data
# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "PerformanceAnalytics",
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
                                                TRUE ~ "Y"))) %>%
  dplyr::filter(substrate %in% c("dimethyl", "TMA"))

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
  summarise_at(vars(activity), mean) %>%
  dplyr::rename(avg_sub = activity) %>%
  arrange(desc(avg_sub)) 
avgsub

combdat <- comb %>%
  right_join(avgsub, ., by = "substrate") %>%
  left_join(., avgdat, by = "org") %>%
  add_count(org) %>%
  dplyr::filter(to_keep == "Y") %>%
  dplyr::select(org, area_sa, avg_activity, volume_sa) %>%
  dplyr::mutate(spec = paste0(substr(org, 1, 1), ". ", word(org, sep = " ", -1))) %>%
  # dplyr::mutate(num_sub = paste0(rep(1:14, each = n), substrate)) %>%
  distinct() 
combdat


my.formula = y ~ x
pdf("output/machine_learning/solvent_pocket_volume_surface_area_plot_dimethyl.pdf", width = 7, height = 6)
pl2 <- ggplot(data = combdat, aes(x = area_sa, y = avg_activity)) +
  geom_smooth(data = combdat, aes(x = area_sa, y = avg_activity), method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  geom_point(data = combdat, aes(x = area_sa, y = avg_activity, size = volume_sa),
             color = "#008080", 
             alpha = ifelse(grepl("atraur|Halobacter|paraconglomeratum", combdat$org), 1, 0.3)) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "none") +
  # scale_color_viridis(option = "plasma") +
  # ggrepel::geom_text_repel(data = subset(combdat, avg_activity > 1.85 | area_sa > 580 | area_sa < 275), segment.size = 0.2,
  #                          # nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
  #                          segment.color = "grey50", #direction = "x",
  #                          aes(label = spec), fontface = "italic") +
  ggrepel::geom_text_repel(data = subset(combdat, grepl("Kytococcus|campestris|combesii|avermitilis", org)),
                           segment.size = 0.2, nudge_x = 2, #force = 50, 
                           segment.color = "gray40", size = 4,
                           # nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
                           color = "gray40",
                           #direction = "x",
                           aes(label = spec), fontface = "italic") +
  ggrepel::geom_text_repel(data = subset(combdat, grepl("atraur|Halobacter|paraconglomeratum", org)),
                           segment.size = 0.2, size = 4, nudge_x = 2, #force = 50,
                           segment.color = "black", 
                           # nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
                           color = "black",
                           #direction = "x",
                           aes(label = spec), fontface = "italic") +
  ylab("Avg. enzyme activity (bulky substrates)") +
  xlab("Solvent-accessible surface area") +
  xlim(c(200, 700)) +
  scale_y_continuous(breaks = seq(0.8, 1.8, by = 0.2)) 
pl2
dev.off()


# Now merge everything...
comb2 <- activity %>%
  dplyr::left_join(., molec_fts, by = "substrate") %>%
  dplyr::left_join(., cavity_fts, by = "org") %>%
  dplyr::left_join(., chem_descr, by = "substrate") %>%
  dplyr::mutate(id = paste0(org, "_", substrate)) %>%
  dplyr::mutate(is_active = as.factor(case_when(activity == 0 ~ "N",
                                                TRUE ~ "Y"))) 

comb2$area_sa <- as.numeric(comb2$area_sa)
comb2$volume_sa <- as.numeric(comb2$volume_sa)

avgdat <- comb2 %>%
 # dplyr::filter(activity != 0) %>%
  group_by(org) %>%
  summarise_at(vars(activity), mean) %>%
  dplyr::rename(avg_activity = activity)

avgsub <- comb2 %>%
  #dplyr::filter(activity != 0) %>%
  group_by(substrate) %>%
  summarise_at(vars(activity), median) %>%
  dplyr::rename(avg_sub = activity) %>%
  arrange(desc(avg_sub)) 
avgsub

combdat2 <- comb2 %>%
  right_join(avgsub, ., by = "substrate") %>%
  left_join(., avgdat, by = "org") %>%
  # dplyr::filter(activity != 0) %>%
  add_count(org) %>%
  # add_count(substrate) %>%
  dplyr::filter(to_keep == "Y") %>%
  #dplyr::select(org, substrate, avg_sub, volume_sa, area_sa, contains("PC"), HybRatio, VABC, nAtomLAC, MLogP, TopoPSA, O, MW, nRotB, avg_activity, activity, n) %>%
  # dplyr::select(org, volume_sa, area_sa, avg_activity, n) %>%
  #dplyr::rename(`Number of pNP substrates accepted` = n) %>%
  # dplyr::rename(`Average activity` = avg_activity) %>%
  dplyr::select(org, area_sa, avg_activity, volume_sa, n) %>%
  dplyr::mutate(spec = paste0(substr(org, 1, 1), ". ", word(org, sep = " ", -1))) %>%
  # dplyr::mutate(num_sub = paste0(rep(1:14, each = n), substrate)) %>%
  distinct()
combdat2


combdat2$avg_activity
combdat2$area_sa

corrdf <- combdat2 %>%
  dplyr::select(-spec, -org)

chart.Correlation(corrdf)


cor.test(corrdf$area_sa, corrdf$avg_activity, method = "pearson")
cor(corrdf$area_sa, corrdf$avg_activity, method = "pearson")

my.formula = y ~ x
set.seed(1234)
pdf("output/machine_learning/solvent_pocket_volume_surface_area_plot_all_substrates_compared.pdf", width = 7, height = 6)
pl1 <- ggplot(data = combdat2, aes(x = area_sa, y = avg_activity)) +
  geom_smooth(data = combdat2, aes(x = area_sa, y = avg_activity), method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  geom_point(data = combdat2, aes(x = area_sa, y = avg_activity, size = volume_sa),
             color = "gray30", 
             alpha = ifelse(grepl("atraur|Halobacter|paraconglomeratum", combdat2$org), 1, 0.3)) +
  theme_pubr(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        legend.title = element_blank()) +
 # scale_color_viridis(option = "plasma") +
  # ggrepel::geom_text_repel(data = subset(combdat, avg_activity > 1.85 | area_sa > 580 | area_sa < 275), segment.size = 0.2,
  #                          # nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
  #                          segment.color = "grey50", #direction = "x",
  #                          aes(label = spec), fontface = "italic") +
  ggrepel::geom_text_repel(data = subset(combdat2, grepl("Kytococcus|campestris|combesii|avermitilis", org)),
                           segment.size = 0.2, size = 4, #force = 30, 
                           segment.color = "gray40", nudge_x = 2,
                           # nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
                           color = "gray40",
                           #direction = "x",
                           aes(label = spec), fontface = "italic") +
  ggrepel::geom_text_repel(data = subset(combdat2, grepl("atraur|Halobacter|paraconglomeratum", org)),
                           segment.size = 0.2, size = 4, #force = 30,
                           segment.color = "black", nudge_x = 2,
                           # nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
                           color = "black",
                           #direction = "x",
                           aes(label = spec), fontface = "italic") +
  ylab("Avg. enzyme activity (all substrates)") +
  xlab("Solvent-accessible surface area") +
  # stat_poly_eq(formula = my.formula,
  #            aes(label = paste( ..rr.label.., sep = "~~~")),
  #            parse = TRUE)
  scale_y_continuous(breaks = seq(0, 2.2, by = 0.2)) +
  xlim(c(200, 700))
pl1
dev.off()


pdf("output/machine_learning/vert_stacked_volume_plot.pdf", height = 9, width = 6)
plot_grid(pl1, pl2, ncol = 1, rel_heights = c(5, 3))
dev.off()

pdf("output/machine_learning/solvent_pocket_volume_surface_area_plot_all_substrates_compared.pdf", width = 7, height = 6)
ggplot(data = combdat2, aes(x = area_sa, y = avg_activity)) +
  geom_smooth(data = combdat2, aes(x = area_sa, y = avg_activity), method = "lm", se = FALSE, color="gray90", formula = my.formula, linetype = 2) +
  geom_point(data = combdat2, aes(x = area_sa, y = avg_activity, size = volume_sa),
             color = "gray30", 
             alpha = ifelse(grepl("atraur|Halobacter|paraconglomeratum", combdat2$org), 1, 0.15)) +
  theme_pubr(base_size = 16) +
  theme(legend.title = element_blank()) +
  # scale_color_viridis(option = "plasma") +
  # ggrepel::geom_text_repel(data = subset(combdat, avg_activity > 1.85 | area_sa > 580 | area_sa < 275), segment.size = 0.2,
  #                          # nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
  #                          segment.color = "grey50", #direction = "x",
  #                          aes(label = spec), fontface = "italic") +
  ggrepel::geom_text_repel(data = subset(combdat2, grepl("Kytococcus|campestris|combesii|avermitilis", org)),
                           segment.size = 0.2, size = 4,
                           segment.color = "gray40", 
                           # nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
                           color = "gray40",
                           #direction = "x",
                           aes(label = spec), fontface = "italic") +
  ggrepel::geom_text_repel(data = subset(combdat2, grepl("atraur|Halobacter|paraconglomeratum", org)),
                           segment.size = 0.2, size = 4,
                           segment.color = "black", 
                           # nudge_y = 20, #subset(combdat, avg_activity > 1.25)$volume_sa,
                           color = "black",
                           #direction = "x",
                           aes(label = spec), fontface = "italic") +
  ylab("Avg. enzyme activity (all substrates)") +
  xlab("Solvent-accessible surface area") +
  # stat_poly_eq(formula = my.formula,
  #            aes(label = paste( ..rr.label.., sep = "~~~")),
  #            parse = TRUE)
  scale_y_continuous(breaks = seq(0, 2.2, by = 0.2)) +
  xlim(c(200, 700))
dev.off()

