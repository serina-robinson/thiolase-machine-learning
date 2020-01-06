# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", "randomcoloR", "ggpubr", "anomalize",
               "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales", "dplyr", "forcats")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

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
  dplyr::select(org, substrate, sub_num, avg_sub, volume_sa, area_sa, contains("PC"), AROMATIC, VABC, nAtomLAC, MLogP, TopoPSA, O, MW, nRotB, avg_activity, activity, n) %>%
  # dplyr::select(org, volume_sa, area_sa, avg_activity, n) %>%
  #dplyr::rename(`Number of pNP substrates accepted` = n) %>%
  #dplyr::rename(`Average activity` = avg_activity) %>%
  group_by(substrate) %>%
  # dplyr::mutate(num_sub = paste0(rep(1:14, each = n), substrate)) %>%
  distinct() 

# Calculate outliers 
outl <- combdat %>%
  group_by(substrate) %>%
  mutate(which_max = max(activity)) %>%
  mutate(is_max = activity == which_max) %>%
  mutate(first_quartile = summary(activity)[2]) %>% # 1st quartile
  mutate(third_quartile = summary(activity)[5]) %>% # 3rd quartile
  mutate(iqr = third_quartile - first_quartile) %>%
  mutate(outlier_thresh = iqr * 1.5) %>%
  mutate(is_outlier = case_when(activity > (third_quartile + outlier_thresh) ~ TRUE,
                                activity < (first_quartile - outlier_thresh) ~ TRUE,
                                TRUE ~ FALSE)) %>%
  dplyr::rename(`Molecular topology index (PC2)` = PC2)

outl$is_max

pal2 <- inferno(14)
outl$PC2
unique(pal2[as.factor(outl$PC2)])



pdf("output/substrate_comparisons/20200105_boxplots_with_means_colored_PC2.pdf", width = 10, height = 6)
ggplot(aes(x = sub_num, y = activity), data = outl) + #fill = "gray80") +
  #aes(x = forcats::fct_reorder(cmpnd, log_slope, fun = mean, .desc = TRUE), y = log_slope), data = maprdat_long) +
  geom_boxplot(aes(fill = `Molecular topology index (PC2)`), outlier.color = NA, colour = "gray60") +
  
  # geom_jitter(data = subset(maprdat_long, is_outlier),
  #             color = "blue", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  #  geom_abline(intercept = 0.5, slope = 0, color = "red") +
  geom_jitter(color = "gray40", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  geom_point(aes(x = sub_num, y = avg_sub),
             fill = "gray40", color = "black", shape = 23, size = 2, alpha = 0.5) +
  
  # geom_segment(aes(x = cmpnd - 1, y = mean_slope, xend = cmpnd, yend = mean_slope), color = "pink") +
  labs_pubr() +
  theme_pubr() +
  ylab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  xlab("Substrate") +
  # ylim(0, 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8)) +
  #coord_flip() +
  scale_fill_viridis(option = "magma") +
  ggrepel::geom_text_repel(data = subset(outl, outl$is_max == TRUE & outl$is_outlier == TRUE), # (combdat, outl$activity > (outl$third_quartile + outl$outlier_thresh)), #segment.size = 0.2,
                          #nudge_y = 1, 
                          nudge_x = 0.6,
                         # nudge_x = 1,#subset(combdat, avg_activity > 1.25)$volume_sa,
                           segment.color = "black", #direction = "x",
                           aes(label = org), fontface = "italic") +
  coord_flip() 
  # theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
  #   axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

dev.off()

pdf("output/substrate_comparisons/20200105_boxplots_with_means_colored_PC5.pdf", width = 10, height = 6)
ggplot(aes(x = sub_num, y = activity), data = outl) + #fill = "gray80") +
  #aes(x = forcats::fct_reorder(cmpnd, log_slope, fun = mean, .desc = TRUE), y = log_slope), data = maprdat_long) +
  geom_boxplot(aes(fill = PC5), outlier.color = NA, colour = "gray60") +
  
  # geom_jitter(data = subset(maprdat_long, is_outlier),
  #             color = "blue", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  #  geom_abline(intercept = 0.5, slope = 0, color = "red") +
  geom_jitter(color = "gray40", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  geom_point(aes(x = sub_num, y = avg_sub),
             fill = "gray40", color = "black", shape = 23, size = 2, alpha = 0.5) +
  
  # geom_segment(aes(x = cmpnd - 1, y = mean_slope, xend = cmpnd, yend = mean_slope), color = "pink") +
  labs_pubr() +
  theme_pubr() +
  ylab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  xlab("Substrate") +
  # ylim(0, 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8)) +
  #coord_flip() +
  scale_fill_viridis(option = "magma") +
  ggrepel::geom_text_repel(data = subset(outl, outl$is_max == TRUE & outl$is_outlier == TRUE), # (combdat, outl$activity > (outl$third_quartile + outl$outlier_thresh)), #segment.size = 0.2,
                           #nudge_y = 1, 
                           nudge_x = 0.6,
                           # nudge_x = 1,#subset(combdat, avg_activity > 1.25)$volume_sa,
                           segment.color = "black", #direction = "x",
                           aes(label = org), fontface = "italic") +
  coord_flip() 
# theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
#   axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

dev.off()


pdf("output/substrate_comparisons/20200105_boxplots_with_means_colored_nRotB.pdf", width = 10, height = 6)
ggplot(aes(x = sub_num, y = activity), data = outl) + #fill = "gray80") +
  #aes(x = forcats::fct_reorder(cmpnd, log_slope, fun = mean, .desc = TRUE), y = log_slope), data = maprdat_long) +
  geom_boxplot(aes(fill = nRotB), outlier.color = NA, colour = "gray60") +
  
  # geom_jitter(data = subset(maprdat_long, is_outlier),
  #             color = "blue", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  #  geom_abline(intercept = 0.5, slope = 0, color = "red") +
  geom_jitter(color = "gray40", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  geom_point(aes(x = sub_num, y = avg_sub),
             fill = "gray40", color = "black", shape = 23, size = 2, alpha = 0.5) +
  
  # geom_segment(aes(x = cmpnd - 1, y = mean_slope, xend = cmpnd, yend = mean_slope), color = "pink") +
  labs_pubr() +
  theme_pubr() +
  ylab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  xlab("Substrate") +
  # ylim(0, 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8)) +
  #coord_flip() +
  scale_fill_viridis(option = "magma") +
  ggrepel::geom_text_repel(data = subset(outl, outl$is_max == TRUE & outl$is_outlier == TRUE), # (combdat, outl$activity > (outl$third_quartile + outl$outlier_thresh)), #segment.size = 0.2,
                           #nudge_y = 1, 
                           nudge_x = 0.6,
                           # nudge_x = 1,#subset(combdat, avg_activity > 1.25)$volume_sa,
                           segment.color = "black", #direction = "x",
                           aes(label = org), fontface = "italic") +
  coord_flip() 
# theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
#   axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

dev.off()

pdf("output/substrate_comparisons/20200105_boxplots_with_means_colored_AROMATIC.pdf", width = 10, height = 6)
ggplot(aes(x = sub_num, y = activity), data = outl) + #fill = "gray80") +
  #aes(x = forcats::fct_reorder(cmpnd, log_slope, fun = mean, .desc = TRUE), y = log_slope), data = maprdat_long) +
  geom_boxplot(aes(fill = AROMATIC), outlier.color = NA, colour = "gray60") +
  
  # geom_jitter(data = subset(maprdat_long, is_outlier),
  #             color = "blue", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  #  geom_abline(intercept = 0.5, slope = 0, color = "red") +
  geom_jitter(color = "gray40", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  geom_point(aes(x = sub_num, y = avg_sub),
             fill = "gray40", color = "black", shape = 23, size = 2, alpha = 0.5) +
  
  # geom_segment(aes(x = cmpnd - 1, y = mean_slope, xend = cmpnd, yend = mean_slope), color = "pink") +
  labs_pubr() +
  theme_pubr() +
  ylab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  xlab("Substrate") +
  # ylim(0, 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8)) +
  #coord_flip() +
  scale_fill_viridis(option = "magma") +
  ggrepel::geom_text_repel(data = subset(outl, outl$is_max == TRUE & outl$is_outlier == TRUE), # (combdat, outl$activity > (outl$third_quartile + outl$outlier_thresh)), #segment.size = 0.2,
                           #nudge_y = 1, 
                           nudge_x = 0.6,
                           # nudge_x = 1,#subset(combdat, avg_activity > 1.25)$volume_sa,
                           segment.color = "black", #direction = "x",
                           aes(label = org), fontface = "italic") +
  coord_flip() 
# theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
#   axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

dev.off()



pdf("output/substrate_comparisons/20200105_boxplots_with_means_sorted_by_cmpnd.pdf", width = 10, height = 6)
ggplot(aes(x = sub_num, y = activity), data = outl) + #fill = "gray80") +
  #aes(x = forcats::fct_reorder(cmpnd, log_slope, fun = mean, .desc = TRUE), y = log_slope), data = maprdat_long) +
  geom_boxplot(outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +
  
  # geom_jitter(data = subset(maprdat_long, is_outlier),
  #             color = "blue", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  #  geom_abline(intercept = 0.5, slope = 0, color = "red") +
  geom_jitter(data = subset(outl, activity > avg_sub),
              color = "firebrick", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  geom_jitter(data = subset(outl, activity < avg_sub),
              color = "blue", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  geom_point(aes(x = sub_num, y = avg_sub),
             fill = "firebrick", color = "black", shape = 23, size = 2, alpha = 0.5) +
  
  # geom_segment(aes(x = cmpnd - 1, y = mean_slope, xend = cmpnd, yend = mean_slope), color = "pink") +
  labs_pubr() +
  theme_pubr() +
  ylab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  xlab("Substrate") +
  # ylim(0, 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8))
        # axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        # axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
dev.off()
