# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", "randomcoloR", "ggpubr", "anomalize",
               "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales", "dplyr", "forcats")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "calculated_slopes.csv"


ll <- fils[grepl(suffix, fils)]
ll <- ll[!grepl("BocPhe|furf|scratch|only|averaged|benzoate|round|2019-11-15_7Ph", ll)] # rep1|rep2|reps|round|
ll
# Remove the 7Ph that was from the bad run (non-normalized)

rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(cmpnd =  paste0(word(word(filename, 2, sep = "output\\/\\/"), 1, sep = "\\/"), "_",
    word(word(filename, 2, sep = "_"), 1, sep = "_")))
                  #word(word(sep = "output//", 2, filename), 1, sep = "\\/"))
  
rawdat$filename[grep("dodecanoate", rawdat$cmpnd)]

# dplyr::mutate(cmpnd = case_when(grepl("C6", filename) ~ "hexanoate",
  #                                 grepl("heptynoate", filename) ~ "heptynoate", 
  #                                 TRUE ~ word(word(filename, 2, sep = "_"), 1, sep = "_")))
# Combine dodecanoate_1, _2, and _3


avged_ctrls <- rawdat %>%
 dplyr::filter(grepl("C6", cmpnd))
avged_ctrls

newdat <- rawdat %>%
  dplyr::filter(!grepl("C6", cmpnd))

# Write to file
maprdat_log <- newdat %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope)) 
maprdat_log

# Combine with c6
maprdat_comb <- maprdat_log %>%
  bind_rows(., avged_ctrls) %>%
  dplyr::select(cmpnd, org, log_slope)
unique(maprdat_comb$cmpnd)

# Calculate outliers 
outl <- maprdat_comb %>%
  group_by(cmpnd) %>%
  mutate(first_quartile = summary(log_slope)[2]) %>% # 1st quartile
  mutate(third_quartile = summary(log_slope)[5]) %>% # 3rd quartile
  mutate(iqr = third_quartile - first_quartile) %>%
  mutate(outlier_thresh = iqr * 1.5) %>%
  mutate(is_outlier = case_when(log_slope > (third_quartile + outlier_thresh) ~ TRUE,
                                log_slope < (first_quartile - outlier_thresh) ~ TRUE,
                                TRUE ~ FALSE))
outl

# Find the average
maprdat_avg <- outl %>%
  #dplyr::filter(is_outlier == FALSE) %>% # remove outliers in calculations
  dplyr::group_by(cmpnd) %>% # this would necessarily be (org, cmpnd) once we have biological replicates!!!!
  dplyr::summarise_each(funs(mean), log_slope) %>%
  dplyr::rename(mean_slope = log_slope)
maprdat_avg

# Combine
maprdat_long <- outl %>%
  left_join(., maprdat_avg, by = "cmpnd")
table(maprdat_long$cmpnd)

maprdat_long$cmpnd <- gsub("_", " ", maprdat_long$cmpnd) 
maprdat_long$cmpnd[grep("C6", maprdat_long$cmpnd)] <- c("hexanoate averaged")

# alph <- word(maprdat_long$cmpnd, sep = " ", 2)
# testr <- maprdat_long[order(alph),]
# maprdat_long$cmpnd <- as.factor(maprdat_long$cmpnd)
# levels(maprdat_long$cmpnd) <- unique(testr$cmpnd)


pdf("output/scratch_output/20191218_boxplots_with_means_sorted_by_date.pdf", width = 10, height = 6)
ggplot(aes(x = cmpnd, y = log_slope), data = maprdat_long) + #fill = "gray80") +
  #aes(x = forcats::fct_reorder(cmpnd, log_slope, fun = mean, .desc = TRUE), y = log_slope), data = maprdat_long) +
  geom_boxplot(outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +

  # geom_jitter(data = subset(maprdat_long, is_outlier),
  #             color = "blue", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
 #  geom_abline(intercept = 0.5, slope = 0, color = "red") +
  geom_jitter(data = subset(maprdat_long, log_slope > mean_slope),
              color = "firebrick", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  geom_jitter(data = subset(maprdat_long, log_slope < mean_slope),
              color = "blue", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  geom_point(aes(x = cmpnd, y = mean_slope),
             fill = "firebrick", color = "black", shape = 23, size = 3) +

  # geom_segment(aes(x = cmpnd - 1, y = mean_slope, xend = cmpnd, yend = mean_slope), color = "pink") +
  labs_pubr() +
  theme_pubr() +
  ylab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  xlab("Substrate") +
 # ylim(0, 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
dev.off()

alph <- word(maprdat_long$cmpnd, sep = " ", 2)
testr <- maprdat_long[order(alph),]
maprdat_long$cmpnd <- as.factor(maprdat_long$cmpnd)
maprdat_long$cmpnd <- factor(maprdat_long$cmpnd,
                       levels = unique(testr$cmpnd), ordered = TRUE)

maprdat_long$cmpnd

pdf("output/scratch_output/20191218_boxplots_with_means_sorted_by_cmpnd.pdf", width = 10, height = 6)
ggplot(aes(x = cmpnd, y = log_slope), data = maprdat_long) + #fill = "gray80") +
  #aes(x = forcats::fct_reorder(cmpnd, log_slope, fun = mean, .desc = TRUE), y = log_slope), data = maprdat_long) +
  geom_boxplot(outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +
  
  # geom_jitter(data = subset(maprdat_long, is_outlier),
  #             color = "blue", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  #  geom_abline(intercept = 0.5, slope = 0, color = "red") +
  geom_jitter(data = subset(maprdat_long, log_slope > mean_slope),
              color = "firebrick", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  geom_jitter(data = subset(maprdat_long, log_slope < mean_slope),
              color = "blue", alpha = 0.5, position=position_jitter(width=.1, height=0)) +
  geom_point(aes(x = cmpnd, y = mean_slope),
             fill = "firebrick", color = "black", shape = 23, size = 2, alpha = 0.5) +
  
  # geom_segment(aes(x = cmpnd - 1, y = mean_slope, xend = cmpnd, yend = mean_slope), color = "pink") +
  labs_pubr() +
  theme_pubr() +
  ylab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  xlab("Substrate") +
  # ylim(0, 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
dev.off()

# hex_which <- maprdat_long %>%
#   dplyr::filter(cmpnd == "hexanoate") %>%
#   # dplyr::filter(log_slope >= mean_slope) %>%
#   dplyr::arrange(desc(log_slope)) %>%
#   dplyr::mutate(log_slope = round(log_slope, 2)) %>%
#   dplyr::mutate(mean_slope = round(mean_slope, 2)) %>%
#   dplyr::select(cmpnd, org, log_slope, mean_slope) %>%
#   write_csv(., "output/C6_control/hexanoate_final_log_slope_rates_for_mBIO_paper.csv")
# hex_which

hept_which <- maprdat_long %>%
  dplyr::filter(cmpnd == "heptynoate") %>%
  dplyr::filter(log_slope >= mean_slope)
hept_which

### Calculate the actives
# maprdat_long$is_active <- ifelse(maprdat_long$log_slope >= 1.5, T, F)
# maprdat_true <- maprdat_count %>%
#   dplyr::filter(is_active == T)
# table(maprdat_true$cmpnd)
# table(maprdat_count$cmpnd)



## Actives only
colnames(maprdat_long)
maprdat_count <- maprdat_long %>%
  # dplyr::filter(is_active) %>%
  group_by(cmpnd) %>%
  add_count() %>%
  dplyr::select(org, cmpnd, log_slope, n)
maprdat_count

# maprdat_uniq <- maprdat_count %>%
#   dplyr::arrange(desc(log_slope)) %>%
#   dplyr::group_by(cmpnd) %>%
#   slice(1)


# Need to match up with master wells 
key <- read_csv("data/72_OleA_masterwell_org_key_updated_for_Megan.csv") 
key$orgs <- gsub("_", " ", key$orgs)
key$orgs[grep("Lysobacter", key$orgs)] <- "Luteimonas tolerans"

maprdat_to_plate <- maprdat_count %>%
  dplyr::filter(grepl("dimethyl|biotin|decanoate|ClPh|Biotin|TMA", cmpnd)) %>%
  inner_join(., key, by = c("org" = "orgs")) %>%
  dplyr::select(org, cmpnd, log_slope, n, master_well)
write_csv(maprdat_to_plate, "output/scratch_output/which_to_plate_for_cmpnds.csv")

