# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", 
               "randomcoloR", "RColorBrewer", "ggplot2", 
               "pheatmap", "viridis", "scales", "car")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "_calculated_slopes.csv"
head(fils)

ll <- fils[grepl(suffix, fils)]
ll <- ll[!grepl("BocPhe|C6|furf|rep1|rep2", ll)]
ll

rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(cmpnd = case_when(grepl("C6", filename) ~ "hexanoate",
                                  TRUE ~ word(word(filename, 2, sep = "_"), 1, sep = "_")))


hist(#breaks = seq(min(rawdat$max_slope), max(rawdat$max_slope), by = 0.1),
     rawdat$max_slope, breaks = 50)

# Untransformed
d <- density(rawdat$max_slope)
plot(d, type="n")
polygon(d, col="red", 
        border = "gray")
qqPlot(rawdat$max_slope)

# Square root
slope_sqrt = sqrt(rawdat$max_slope)
hist(slope_sqrt, breaks = 50)
qqPlot(slope_sqrt)

# Cube root
slope_cub = sign(rawdat$max_slope) * abs(rawdat$max_slope)^(1/3)
hist(slope_cub)
qqPlot(slope_cub)

# Natural log
slope_log = log(rawdat$max_slope)
hist(slope_log)
qqPlot(slope_log)

# Log base 10
slope_log10 = log10(rawdat$max_slope)
hist(slope_log10)
qqPlot(slope_log10)
