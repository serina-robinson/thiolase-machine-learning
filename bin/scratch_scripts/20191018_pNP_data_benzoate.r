# Install packages
pacman::p_load("tidyverse", "readxl")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the template
# Read in the template
temp <- read_excel("data/Plate Set up SynBio Paper.xlsx", skip = 3) %>%
  janitor::clean_names() %>%
  dplyr::select(-x1, -x11, -x12) %>%
  as.matrix() %>%
  t %>%
  as.vector()
temp

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
head(dat)
head(dat[,2:ncol(dat)])
which.max(colSums(sapply(dat[,2:ncol(dat)], as.numeric))) # E4


colnames(dat) <- c("time", temp)
colnames(dat)
dat2[which.max(dat2$value),]

oleA <- read_csv("data/72_OleA_masterwell_org_key.csv")
newnams <- left_join(data.frame(temp, stringsAsFactors = F), oleA, by = c("temp" = "master_well")) %>% 
  dplyr::mutate(orgs = case_when(is.na(orgs) ~ temp,
                                 TRUE ~ orgs)) %>%
  dplyr::select(temp, orgs)

newnams$orgs[grep("E4", newnams$temp)]


colnames(dat) <- c("time", newnams$orgs)


dat2 <- dat %>%
  reshape2::melt(., id = 'time') %>%
  dplyr::mutate(value = as.numeric(value)) %>%
#  dplyr::mutate(winners = case_when(grepl("e4", variable) ~ as.character(variable),
#                                                                       TRUE ~ " inactive")) %>%
  dplyr::mutate(winners = case_when(grepl("pNP|Pet|XC|Kocuria varians", variable) ~ as.character(variable),
                                    TRUE ~ " inactive"))
  
dat2$value
table(dat2$winners)

# Plotting 
#pal2 <- distinctColorPalette(length(unique(dat2$winners)))
pal <- colorRampPalette(brewer.pal(8,"Set1"))(8)
pal2 <- c("gray80", pal[c(1:5, 8)], "tan4", "salmon", "gold1", "black")

pdf("output/20191019_benzoate_screening_assay.pdf", width = 20, height = 10)
pl <- ggplot(dat2, aes(x=time, y=value, color = winners)) +  #color=variable)) +
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
        legend.position="top") +
  scale_color_manual(values=pal2)
pl
dev.off()




