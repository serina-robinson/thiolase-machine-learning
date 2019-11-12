# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", "randomcoloR", 
               "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "calculated_slopes"
fils

ll <- fils[grepl(suffix, fils)]
scratch <- read_csv(file.path(ll[grepl("heptynoate_average", ll)]))
scratch
ll <- ll[!grepl("BocPhe|C6|furf|rep1|rep2|scratch|reps", ll)]
# ll2 <- c(scratch, ll)
# ll2

rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(cmpnd = case_when(grepl("C6", filename) ~ "hexanoate",
                                  TRUE ~ word(word(filename, 2, sep = "_"), 1, sep = "_")))

# Write to file
maprdat_long <- rawdat %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope)) %>%
  bind_rows(scratch)

# Find the average
maprdat_avg <- maprdat_long %>%
  dplyr::group_by(cmpnd) %>% # this would necessarily be (org, cmpnd) once we have biological replicates!!!!
  dplyr::summarise_each(funs(mean), log_slope)

maprdat_long$is_active <- ifelse(maprdat_long$log_slope >= 1.5, T, F)


## Actives only
maprdat_count <- maprdat_long %>%
  # dplyr::filter(is_active) %>%
  group_by(cmpnd, is_active) %>%
  add_count() %>%
  dplyr::select(org, cmpnd, log_slope, is_active, n)

maprdat_true <- maprdat_count %>%
  dplyr::filter(is_active == T)
table(maprdat_true$cmpnd)

maprdat_uniq <- maprdat_count %>%
  dplyr::filter(cmpnd %in% c("heptynoate", "Butoxy", "Azido")) %>%
  dplyr::arrange(desc(log_slope)) %>%
  dplyr::group_by(cmpnd, is_active) %>%
  dplyr::slice(1)

maprdat_uniq
sum(maprdat_uniq$n)

# There are 245 positive hits, which is less than 3 plates
sum(maprdat_uniq$n[maprdat_uniq$is_active == "TRUE"])
# There are 181 'positive' hits using our cut-off
# That is less than two plates!!

pdf("output/scratch_output/boxplots_with_avgs.pdf", width = 10)
ggplot(aes(x = cmpnd, y = log_slope), data = maprdat_long) +
  # geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = subset(maprdat_long, is_active), color = "blue",
              position=position_jitter(width=.1, height=0)) +
  geom_jitter(data = subset(maprdat_long, !is_active), color = "red",
              position=position_jitter(width=.1, height=0)) +
 # geom_jitter(position=position_jitter(width=.1, height=0)) +
  geom_abline(intercept = 1.5, slope = 0, color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))
dev.off()
