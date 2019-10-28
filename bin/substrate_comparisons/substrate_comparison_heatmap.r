# Pipeline for files with # Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer", "ggplot2", "pheatmap", "viridis", "scales")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("output/", recursive=TRUE, full.names = T)
suffix <- "_slopes_calculated.csv"


ll <- fils[grepl(suffix, fils)]
rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.)


maprdat_long  <- rawdat %>%
  dplyr::filter(!grepl("Pet28", variable))
head(maprdat_long)

maprdat_wide <- dcast(maprdat_long, variable ~ cmpnd, value.var = "minutes") 
maprdat_wide

maprdat_mat <- maprdat_wide %>%
  dplyr::select(-variable) %>%
  as.matrix()

summary(maprdat_long$minutes)

maprdat_mat[is.na(maprdat_mat)]

# sum(is.infinite(maprdat_mat))
maprdat_mat[(maprdat_mat < 0)] <- 0
maprdat_mat

rownames(maprdat_mat) <- maprdat_wide$variable
maprdat_mat <- maprdat_mat #%>%
#  t()

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(maprdat_mat, n = 11)
mat_breaks
labs <- as.character(scientific(mat_breaks, digit = 3))

pdf("output/substrate_comparison_heatmap_clustered.pdf") #width = 10, height = 8)
pheatmap(
         legend_labels = labs,
         maprdat_mat, 
         breaks = mat_breaks,
         color = inferno(5),
         annotation_names_row = T, 
         fontsize_col = 10, 
         fontsize_row = 6, 
         annotation_names_col =  T)
dev.off()

# Ggplot2 heatmap
# pdf("output/substrate_comparison_heatmap.pdf")
# ggplot(maprdat_long, aes(cmpnd, variable)) +
#   geom_tile(aes(fill = value)) + 
#  #  geom_text(aes(label = round(value, 1))) +
#   scale_fill_gradient(low = "white", high = "royalblue") +
#   theme(axis.text.x = element_text(angle = 90))
# dev.off()

