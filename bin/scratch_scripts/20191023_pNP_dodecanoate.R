# Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the template
rawtemp <- read_excel("data/Plate Set up SynBio Paper.xlsx", skip = 3)
rawtemp
temp <- read_excel("data/Plate Set up SynBio Paper.xlsx", skip = 3) %>%
  janitor::clean_names() %>%
  dplyr::select(-x1, -x11, -x12) %>%
  as.matrix() %>%
  t %>%
  as.vector()
temp

# Read in the TMA
tmafils <- list.files("data/", pattern = "dodecanoate", full.names = T)
tmafils <- tmafils[!grepl("~", tmafils)]
tmafils

tma1 <- read_excel(tmafils[2]) %>%
  bind_rows(read_excel(tmafils[2])) %>% #  tmafils[2]
  janitor::clean_names() %>%
  dplyr::select(-temperature_c) %>%
  dplyr::select_if(~ !any(is.na(.))) %>%
  dplyr::select_if(!grepl("11|12", colnames(.))) # remove columns 11 & 12

head(tma1)

oleA <- read_csv("data/72_OleA_masterwell_org_key.csv")
head(oleA)
temp
newnams <- left_join(data.frame(temp, stringsAsFactors = F), oleA, by = c("temp" = "master_well")) %>% 
  dplyr::mutate(orgs = case_when(is.na(orgs) ~ temp,
                     TRUE ~ orgs)) %>%
  dplyr::select(temp, orgs)

colnames(tma1) <- c("time", newnams$orgs)
colnames(tma1) 

dat2 <- tma1 %>%
  reshape2::melt(., id = 'time') %>%
  dplyr::mutate(value = as.numeric(value))%>%
  group_by(variable, time) %>%
  summarise_each(funs(mean, sd), value) %>%
  dplyr::filter(!grepl("pNP", variable))

max(dat2$time)

head(dat2)
# Find the winners
dat3 <- dat2 %>%
  dplyr::filter(grepl(max(time), time)) %>%
  dplyr::mutate(winners = case_when(mean >= 0.5 ~ variable)) %>%
  dplyr::mutate(winners = case_when(is.na(winners) ~ "inactive",
                                    TRUE ~ as.character(winners)))
tofind <- dat3$winners[dat3$winners != "inactive"] # 19 winners
tofind

dat4 <- dat2 %>%
  dplyr::mutate(winners = case_when(variable %in% tofind ~ as.character(variable),
                          TRUE ~ " inactive")) 
  #dplyr::filter(variable %in% c("XC", "Pet28"))

# pal2 <- distinctColorPalette(length(unique(dat4$winners)))
#rawpal <- read_csv("data/OleA_palette_key.csv")
#pal <- rawpal$pal2[!rawpal$pal2 %in% c("seagreen1", "cyan3", "gold1", "wheat3", "orchid1", "blueviolet", "#7B6E4F")]
pal <- colorRampPalette(brewer.pal(12,"Paired"))(12)
pal2 <- c("gray80", pal[1:10], "black")
length(pal2)

pdf("output/20191023_dodecanoate_JGI_genes_no_errorbars_rep2.pdf", width = 14, height = 10)
pl <- ggplot(dat4, aes(x=time, y=mean, color=winners)) +
  geom_point() +
  labs(y="Absorbance (410 nm)", x="Time (minutes)") +
  # geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 20),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position="top") + 
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal2) 
pl
dev.off()

# Read in the key of OleAs


