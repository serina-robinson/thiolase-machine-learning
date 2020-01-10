# Install packages
pacman::p_load("tidyverse", "maditr", "readxl", "randomcoloR", 
               "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "calculated_slopes"

ll <- fils[grepl(suffix, fils)]
ll <- ll[!grepl("BocPhe|furf|scratch|only|averaged|benzoate|round|2019-11-15_7Ph", ll)] 
ll # 42 files corresponding to 3 x 15 substrates (but hexanoate is already averaged from Megan's analysis, so only 1 file there,
# then there is 

rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) %>%
  dplyr::mutate(cmpnd =  paste0(word(word(filename, 2, sep = "_"), 1, sep = "_"))) %>%
  dplyr::mutate(cmpnd = case_when(grepl("C6", filename) ~ "hexanoate",
                                  grepl("7Ph", cmpnd) ~ "7Ph heptanoate",
                                  grepl("Azido", cmpnd) ~ "azido",
                                  grepl("Butoxy", cmpnd) ~ "butoxy",
                                  grepl("Biotin", cmpnd) ~ "biotin",
                                  grepl("decanoate rep-1", cmpnd) ~ "decanoate",
                                  grepl("decanoate rep-2", cmpnd) ~ "decanoate",
                                  grepl("decanoate ", cmpnd) ~ "decanoate",
                                  grepl("ClPh", cmpnd) ~ "ClPh propionate",
                                  TRUE ~ cmpnd)) %>%
  group_by(filename) %>%
  dplyr::slice(1) 

table(rawdat$filename, rawdat$cmpnd) %>%
  colSums()
