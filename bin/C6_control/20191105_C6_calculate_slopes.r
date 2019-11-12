# Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer", "ggplot2", "broom", "maditr")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

### Fill in your parameters of interest ###
### 2019-11-05 re-analysis of Megan's C6 data
cmpnd <- "C6"
date <- "2019-11-11"

# Read in Megan's data
set1 <- read_csv("data/C6_control/2019-11-08_whole_cell_JGI_1_melted_bio_rep_2.csv")
unique(set1$variable)
set2 <- read_csv("data/C6_control/2019-09-09_whole_cell_JGI_2_melted_rep_2.csv")
unique(set2$variable)

rawdat <- read_csv("data/C6_control/2019-11-08_whole_cell_JGI_1_melted_bio_rep_2.csv") %>%
   bind_rows(., read_csv("data/C6_control/2019-09-09_whole_cell_JGI_2_melted_rep_2.csv"))
table(rawdat$variable)
 
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
a <- dat5 %>%
  dplyr::mutate(variable = as.character(variable)) %>%
  dplyr::mutate(minutes = as.numeric(time)) %>%
  dplyr::group_by(variable) %>%
  dplyr::slice(0:45) %>%
  dplyr::ungroup() %>%
  dplyr::select(variable, minutes, mean)
head(a)

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

# Calculate slopes for each organisms organism
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
pal2 <- c("gray80", "dodgerblue", "goldenrod", "chartreuse", pal[c(1, 3:5, 8)], "blue", "gold1", "black", distinctColorPalette(60))

pdf(paste0("output/C6_control/", date, "_", cmpnd, "_JGI_genes_test_assessment_round2.pdf"),  width = 16, height = 14)
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

write_csv(slope_merg, paste0("output/C6_control/", date, "_", cmpnd, "_all_data_calculated_slopes_round2.csv"))
