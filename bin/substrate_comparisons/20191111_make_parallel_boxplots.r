# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", "randomcoloR", "ggpubr",
               "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "calculated_slopes.csv"
fils

ll <- fils[grepl(suffix, fils)]
ll
ll <- ll[!grepl("BocPhe|C6|furf|rep1|rep2|scratch|reps", ll)]
ll

rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(cmpnd = case_when(grepl("C6", filename) ~ "hexanoate",
                                  TRUE ~ word(word(filename, 2, sep = "_"), 1, sep = "_")))

# Write to file
maprdat_log <- rawdat %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope))



# Find the average
maprdat_avg <- maprdat_log %>%
  dplyr::group_by(cmpnd) %>% # this would necessarily be (org, cmpnd) once we have biological replicates!!!!
  dplyr::summarise_each(funs(mean), log_slope) %>%
  dplyr::rename(mean_slope = log_slope)
maprdat_avg

# Combine
maprdat_long <- maprdat_log %>%
  left_join(., maprdat_avg, by = "cmpnd")

maprdat_long$is_active <- ifelse(maprdat_long$log_slope >= 1.5, T, F)


## Actives only
maprdat_count <- maprdat_long %>%
  # dplyr::filter(is_active) %>%
  group_by(cmpnd, is_active) %>%
  add_count() %>%
  dplyr::select(filename, org, cmpnd, log_slope, is_active, n)


maprdat_true <- maprdat_count %>%
  dplyr::filter(is_active == T)
table(maprdat_true$cmpnd)
table(maprdat_count$cmpnd)


maprdat_uniq <- maprdat_count %>%
  dplyr::arrange(desc(log_slope)) %>%
  dplyr::group_by(cmpnd, is_active) %>%
  slice(1)
maprdat_uniq

sum(maprdat_uniq$n)
# There are 245 positive hits, which is less than 3 plates
sum(maprdat_uniq$n[maprdat_uniq$is_active == "TRUE"])
# There are 181 'positive' hits using our cut-off
# That is less than two plates!!


maprdat_long$mean_slope

#pdf("output/scratch_output/boxplots_with_means.pdf")
ggplot(aes(x = cmpnd, y = log_slope), data = maprdat_long) +
  # geom_violin() +
  geom_boxplot(geom_boxplot(outlier.colour = "red", outlier.shape = 1)) +
  geom_jitter(data = subset(maprdat_long, is_active), color = "blue",
              position=position_jitter(width=.1, height=0)) +
  geom_jitter(data = subset(maprdat_long, !is_active), color = "red",
              position=position_jitter(width=.1, height=0)) +
 # geom_jitter(position=position_jitter(width=.1, height=0)) +
  geom_abline(intercept = 1.5, slope = 0, color = "red") +
  geom_point(aes(x = cmpnd, y = mean_slope),
             colour="green", shape=18, size=3) +
  theme_pubr()
#dev.off()
