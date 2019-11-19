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
  right_join(., a, by = c("org" = "variable")) #%>%
  #left_join(., a, by = c("org" = "variable")) %>% # to exclude inactive ones
  # dplyr::mutate(winners = case_when(is.na(max_slope) ~ " inactive",
  #                                   TRUE ~ org))


# Set random number seed for random color palette
set.seed(1234)
pal <- colorRampPalette(brewer.pal(8,"Set1"))(8)
pal2 <- c("dodgerblue", "goldenrod", "chartreuse", pal[c(1, 3:5, 8)], "blue", "gold1", "black", "#FF0000", "#0066CC","#008080", 
          "#CC99FF", "#FFCC00", 
          "#666699", "#FFB6C1", "#DB7093", "#333333", "#00FFFF", 
          "#CCCCFF", "#800080", "#99CC00", "#FF6600","#7FFFD4", "#0000FF", "#A52A2A", "#7FFF00", "gray80")

# Read in the active slopes
actives <- read_csv("output/C6_control/hexanoate_final_log_slope_rates_for_mBIO_paper.csv") %>%
  dplyr::filter(log_slope >= mean_slope) %>%
  dplyr::mutate(active = org) %>%
  dplyr::mutate(colr = pal2[1:length(org)])
head(actives)

#halo <-# merg_all[grep(paste0(org, collapse = "|"), merg_all$org),]
keep <- merg_all %>%
  inner_join(., actives, by = "org")
head(keep)

pdf(paste0("output/C6_control/", date, "_", cmpnd, "_JGI_genes_test_assessment_26_active_only_rep2.pdf"),  width = 10, height = 10)
pl <- ggplot(keep,  aes(x = minutes, y = mean)) + 
  geom_point(color = keep$colr) +
  geom_abline(slope = unique(keep$max_slope), intercept = unique(keep$intercept), color = unique(keep$colr)) +
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
  scale_color_manual(values = unique(keep$colr))
pl
dev.off() 

legend_txt <- unique(keep$org)[order(unique(keep$max_slope), decreasing = T)]
color_legend <-  unique(keep$colr)[order(unique(keep$max_slope), decreasing = T)]
legend_txt

pdf(paste0("output/C6_control/", date, "_", cmpnd, "_JGI_genes_test_assessment_26_active_only_rep2_legend.pdf"), width = 15, height = 10)
plot.new()
legend("bottomright", pch=19,
       legend = legend_txt, border = NULL, ncol = 4,
       col = color_legend, bty = "n", text.col="black")
dev.off()

tow <- keep %>%
  group_by(org) %>%
  dplyr::slice(1) %>%
  dplyr::select(org, max_slope, cmpnd) %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope)) %>%
  dplyr::select(org, log_slope) %>%
  dplyr::arrange(desc(log_slope)) %>%
  dplyr::mutate(log_slope = round(log_slope, 2)) 

colnames(tow) <- c("Source organism", "Enzyme activity")
write_csv(tow, paste0("output/C6_control/", date, "_", cmpnd, "_JGI_genes_test_assessment_26_active_only_rep2_slopes.csv"), col_names = T)
# write_csv(slope_merg, paste0("output/C6_control/", date, "_", cmpnd, "_all_data_calculated_slopes_round2.csv"))