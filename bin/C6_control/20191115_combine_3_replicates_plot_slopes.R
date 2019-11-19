# Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer", "ggplot2", "broom", "maditr")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

### Fill in your parameters of interest ###
### 2019-11-05 re-analysis of Megan's C6 data
cmpnd <- "C6"
date <- "2019-11-15"

# Read in Megan's data
rep1_set1 <- read_csv("data/C6_control/2019-09-21_whole_cell_JGI_1_melted.csv")
rep1_set2 <- read_csv("data/C6_control/2019-08-31_whole_cell_JGI_2_melted.csv")
rep2_set1 <- read_csv("data/C6_control/2019-11-08_whole_cell_JGI_1_melted_bio_rep_2.csv")
rep2_set2 <- read_csv("data/C6_control/2019-09-09_whole_cell_JGI_2_melted_rep_2.csv")

rawdat <- rep1_set1 %>%
  bind_rows(rep1_set2, rep2_set1, rep2_set2)
 
oleA <- read_csv("data/72_OleA_masterwell_org_key_updated_for_Megan.csv") %>%
  dplyr::mutate(nam_match = case_when(master_well == "OleA" ~ master_well,
                                      TRUE ~ paste0(word(master_well, sep = "-", 2), ".", word(master_well, sep = "-", 1)))) 

oleA$orgs <- gsub("Lysobacter", "Luteimonas", oleA$orgs)
oleA$orgs <- gsub("_", " ", oleA$orgs)

dat5 <- left_join(data.frame(rawdat, stringsAsFactors = F), oleA, 
                     by = c("variable" = "nam_match")) %>%
  dplyr::filter(variable != "E6.2") %>%
  dplyr::mutate(orgs = case_when(is.na(orgs) ~ variable,
                                 TRUE ~ orgs)) %>%
  dplyr::mutate(oldvar = variable) %>%
  dplyr::mutate(variable = orgs) %>%
  dplyr::rename(time = Time) %>%
  dplyr::mutate(mean = value) %>%
  dplyr::select(time, variable, mean, oldvar)


# Calculate slopes
raw_a <- dat5 %>%
  dplyr::mutate(variable = as.character(variable)) %>%
  dplyr::mutate(minutes = as.numeric(time)) %>%
  dplyr::filter(minutes <= 45) %>%
  # dplyr::group_by(variable) %>%
  # dplyr::slice(0:46) %>%
  # dplyr::ungroup() %>%
  dplyr::select(variable, minutes, mean)
table(raw_a$variable)

# Read in and bind rows 
raw_b <- read_csv("output/C6_control/2019-11-08_whole_cell_hexanoate_melted_bio_rep_3.csv")
raw_c <- raw_a %>%
  bind_rows(raw_b)

table(raw_c$variable, raw_c$minutes)

a <- raw_c %>%
  dplyr::filter(minutes <= 45) %>%
  group_by(variable, minutes) %>%
  dplyr::summarise_each(funs(sd, mean), mean)

## Function that calculates slopes
slopes <- function(d) { 
  m <- lm(mean ~ minutes, as.data.frame(d, stringsAsFactors = F))
  summ <- summary(m)
  r2 <- summ$r.squared
  intercept <- coef(m)[1]
  slope <- coef(m)[2]
  return(list(org = d$variable[1], r2 = r2, slope = slope, intercept = intercept))
}

## The number of frames to take the windowed slope of
windowsize <- 15 # Note had to change window size to 5 for butoxy because Kytococcus was so fast

# Calculate slopes for each organism
orgs <- unique(a$variable)

res <- list()
for(i in 1:length(orgs)) {
  tmp <- a[a$variable == orgs[i],]
  res[[i]] <- do.call(rbind.data.frame,lapply(seq(dim(tmp)[1] - windowsize),
                                              function(x) slopes(tmp[x:(x + windowsize),])))
}

names(res) <- orgs
resl <- plyr::ldply(res, data.frame)
resll <- do.call(rbind.data.frame, res)
resll$org

# Find max slope for each organism
resmax <- resl %>%
  dplyr::filter(r2 >= 0.9) %>% # make sure R^2 is above or equal to 0.9
  group_by(org) %>%
  summarise_each(funs(max_slope = max), slope) %>%
  dplyr::filter(max_slope > 0)
resmax

# Merge with the original dataset
slope_merg <- resmax %>%
  inner_join(., resl, by = "org") %>%
  dplyr::filter(r2 >= 0.9) %>%
  group_by(org) %>%
  dplyr::filter(slope == max(slope)) %>%
  dplyr::select(org, max_slope, r2, intercept)

# Plot winners on graph
merg_all <- slope_merg %>% 
  right_join(., a, by = c("org" = "variable")) %>%
  #left_join(., a, by = c("org" = "variable")) %>% # to exclude inactive ones
  dplyr::mutate(winners = case_when(is.na(max_slope) ~ " inactive",
                                    TRUE ~ org))

pal <- colorRampPalette(brewer.pal(8,"Set1"))(8)
pal2 <- c("gray80", "dodgerblue", "goldenrod", "chartreuse", pal[c(1, 3:5, 8)], "blue", "gold1", "black", distinctColorPalette(70))

pdf(paste0("output/C6_control/", date, "_", cmpnd, "_JGI_genes_all_three_replicates_combined.pdf"),  width = 20, height = 14)
pl <- ggplot(merg_all,  aes(x = minutes, y = mean, color = winners)) + 
  geom_point(alpha = ifelse(merg_all$winners == " inactive", 0.2, 1)) +
  geom_abline(slope = unique(merg_all$max_slope), intercept = unique(merg_all$intercept), color = pal2[1:length(unique(merg_all$max_slope))]) +
  labs(y = "nmol pNP", x = "Time (minutes)") +
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 16),
        legend.position = "top",
        legend.key= element_rect(fill=NA, color=NA)) +
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal2) 
# legend.position = "none")
pl
dev.off() 

# which_max <- merg_all %>%
#   dplyr::filter(max_slope >= 1.63)

# write_csv(slope_merg, paste0("output/C6_control/", date, "_", cmpnd, "_all_data_calculated_slopes_round2.csv"))

# Read in the active ones only and just plot those
actives <- read_csv("output/C6_control/hexanoate_final_log_slope_rates_for_mBIO_paper.csv") %>%
  dplyr::filter(log_slope >= mean_slope) %>%
  dplyr::mutate(active = org)# 26 actives
head(actives)

active_merg <- merg_all %>%
  left_join(., actives, by = "org") %>%
  dplyr::filter(!is.na(active)) %>%
  dplyr::select(org, max_slope, intercept, mean, minutes, cmpnd, active)


pdf(paste0("output/C6_control/", date, "_", cmpnd, "_JGI_genes_all_three_replicates_combined_only_active.pdf"),  width = 16, height = 14)
pl <- ggplot(active_merg,  aes(x = minutes, y = mean, color = active)) + 
  geom_point(#(color = ifelse(is.na(active_merg$active), "gray30", pal2[]), 
             alpha = ifelse(is.na(active_merg$active), 0, 1)) +
  geom_abline(slope = unique(active_merg$max_slope), intercept = unique(active_merg$intercept), 
              color = pal2[1:length(unique(active_merg$max_slope))]) +
  labs(y = "nmol pNP", x = "Time (minutes)") +
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 16),
        legend.position = "top",
        legend.key= element_rect(fill=NA, color=NA)) +
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal2) 
# legend.position = "none")
pl
dev.off() 


calc_slopes <- active_merg %>%
  group_by(org) %>%
  dplyr::slice(1) %>%
  dplyr::select(org, max_slope, cmpnd) %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope)) %>%
  write_csv(., "output/calculated_slopes_averaged.csv")

pdf(paste0("output/C6_control/", date, "_", cmpnd, "_JGI_genes_activity_barplot.pdf"),  width = 16, height = 14)
pl <- ggplot(  aes(x = minutes, y = mean, color = active)) + 
  geom_bar()


