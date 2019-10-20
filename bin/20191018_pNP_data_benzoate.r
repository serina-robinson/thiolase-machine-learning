# Install packages
pacman::p_load("tidyverse", "readxl")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the template
temp <- read_excel("data/Plate Set up SynBio Paper.xlsx", skip = 3) %>%
  janitor::clean_names() %>%
  dplyr::select(-x1, -x11, -x12) %>%
  as.matrix() %>%
  as.vector()
head(temp)

# Read in benzoate
benz <- list.files("data/", pattern = "benz", full.names = T)
benz

# Read in the data
dat <- read_excel(benz) %>%
  janitor::clean_names() %>%
  dplyr::select(-temperature_c) %>%
  dplyr::slice(-62, -63, -64) %>%
  dplyr::select_if(~ !any(is.na(.)))
  # dplyr::select(1:10) %>%

colnames(dat) <- c("time", temp)


dat2 <- dat %>%
  reshape2::melt(., id = 'time') %>%
  dplyr::mutate(value = as.numeric(value))
dat2$value

# Plotting 
pal2 <- distinctColorPalette(length(unique(dat2$variable)))


pdf("output/20191019_benzoate_screening_assay.pdf", width = 10, height = 10)
pl <- ggplot(dat2, aes(x=time, y=value)) +  #color=variable)) +
  geom_point() +
  labs(y="Absorbance (410 nm)", x="Time (minutes)") +
  #geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 20),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position="bottom") +
  scale_color_manual(values=pal2)
pl
dev.off()




