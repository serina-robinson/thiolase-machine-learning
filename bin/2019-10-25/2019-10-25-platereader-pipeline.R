# Pipeline for files with # Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer", "ggplot2", "broom", "maditr")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

### Fill in your parameters of interest ###
### 2019-10-18
# cmpnd <- "dodecanoate"
# cmpnd <- "TMA"
### 2019-10-24
cmpnd <- "7Ph heptanoate"
# cmpnd <- "ClPh propionate"
date <- "2019-10-25"


# Read in the plate template
temp <- read_excel(paste0("data/", date, "/Plate Set up SynBio Paper.xlsx"), skip = 3) %>%
  janitor::clean_names() %>%
  dplyr::select(-x1, -x12) %>% # NOTE THAT X11 WILL BE FILLED 
  as.matrix() %>%
  t %>%
  as.vector()
temp

# Read in the raw data 
tmafils <- list.files(paste0("data/", date, "/"), pattern = cmpnd, full.names = T)
tmafils <- tmafils[!grepl("~", tmafils)] # remove any temporary files
tmafils

## Need to apply to both plates
#### Apply normalize_all to the complete file names
normalize_all <- function(x) { 
  tma1 <- read_excel(x) %>%
    janitor::clean_names() %>%
    dplyr::select(-temperature_c) %>% # remove temp
    #dplyr::select_if(~ !any(is.na(.))) %>% # remove NAs
    dplyr::select_if(!grepl("12", colnames(.))) # remove columns 12 (empty) # TODO change this to keep pet28
  
  # Read in the organism names
  oleA <- read_csv("data/72_OleA_masterwell_org_key.csv")
  oleA$orgs <- gsub("Lysobacter", "Luteimonas", oleA$orgs)
  oleA$orgs <- gsub("_", " ", oleA$orgs)
  newnams <- left_join(data.frame(temp, stringsAsFactors = F), oleA, 
                       by = c("temp" = "master_well")) %>%
    dplyr::mutate(orgs = case_when(is.na(orgs) ~ temp,
                                   TRUE ~ orgs)) %>%
    dplyr::select(temp, orgs)

  newnams$orgs

  # Rename columns to match organisms
  colnames(tma1) <- c("time", newnams$orgs)
  
  # Convert from wide to long format
  dat0 <- tma1 %>%
    reshape2::melt(., id = 'time') %>%
    dplyr::mutate(value = as.numeric(value))
  
  ## TODO: average 3 Pet28s
  pet28 <- dat0 %>%
    dplyr::filter(grepl("Pet28", variable)) %>% # TODO change this to Pet28
    group_by(time) %>%
    summarise_each(funs(mean, sd), value) %>%
    dplyr::pull(mean)
  
  # Remove pNP standards
  dat1 <- dat0 %>%
    dplyr::filter(!grepl("pNP", variable))
  
  pnp_rem <- dat0 %>% 
    dplyr::filter(grepl("pNP", variable))
  
  # Subtract the pET 28b+ empty vector control
  norm_pet <- function(x, na.rm = FALSE) (x - pet28)
  
  dat <- dat1 %>%
    group_by(variable) %>%
    mutate_at(c("value"), norm_pet) %>%
    bind_rows(pnp_rem)

  return(dat)
}

# Normalize all
res <- lapply(tmafils, normalize_all)
head(res[[1]])
head(res[[2]])
table(res[[1]]$value == res[[2]]$value) # check values are different

resbind <- res[[1]] %>%
  bind_rows(res[[2]]) %>%
  dplyr::filter(!grepl("Nothing_", variable))

resbind <- resbind[complete.cases(resbind),] # Removes the "NOTHING's"

# Now normalize to the pNP standard curve
head(pNPs)
pNPs <- resbind %>% 
  dplyr::filter(grepl("pNP", variable)) %>%
  dplyr::mutate(µL = as.numeric(word(variable, sep = " ", 1))) %>%
  dplyr::mutate(mM = µL * (4/200)) %>% # 4 mM stock solution, 200 µL final well volume
  dplyr::mutate(nM = mM * 1e6) %>%
  dplyr::mutate(nmol = nM * (1/1000) * (1/1000) * 200) # nmoles = nmoles/L * 1L/1000 mL * 1ml/1000µL * 200 µL (final volume)
 # dplyr::filter(time == min(time))
pNPs

pNP_fit <- lm(value ~ nmol, data = pNPs)
pNP_fit$coefficients


# Plot the standard curve
pdf(paste0("output/", date, "/", date, "_", cmpnd, "_standard_curve.pdf"))
pl <- ggplot(pNPs,  aes(x = nmol, y = value, color = time)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y="Absorbance (410 nm)", x="nmol pNP") +
  #geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 8),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position="top") + 
  guides(shape = guide_legend(override.aes = list(size = 10)))
pl
dev.off()
pl

# Calculate slope and intercept of pNP standard curve
pNP_fit$coefficients
b <- pNP_fit$coefficients[1]
m <- pNP_fit$coefficients[2]

dat2 <- resbind %>%
  dplyr::filter(!grepl("pNP", variable))
head(dat2)
dat3 <- dat2 %>%
  dplyr::mutate(nmols_pNP = (value - b)/m) %>%
  group_by(variable, time) %>%
  summarise_each(funs(mean, sd), nmols_pNP) 
head(dat3)


# Find the 'winners' i.e. those with a final end point above -0.0005
dat4 <- dat3 %>%
  #dplyr::filter(grepl(max(time), time)) %>%
  dplyr::mutate(winners = case_when(any(mean > 0) ~ variable)) %>%
  dplyr::mutate(winners = case_when(is.na(winners) ~ " inactive",
                                    TRUE ~ as.character(winners)))
dat4$winners
tofind <- unique(dat4$winners)[unique(dat4$winners) != " inactive"] # 9 winners
tofind

dat5 <- dat3 %>%
  dplyr::mutate(winners = case_when(variable %in% tofind ~ as.character(variable),
                                    TRUE ~ " inactive")) 

# pal2 <- distinctColorPalette(length(unique(dat4$winners)))
#rawpal <- read_csv("data/OleA_palette_key.csv")
#pal <- rawpal$pal2[!rawpal$pal2 %in% c("seagreen1", "cyan3", "gold1", "wheat3", "orchid1", "blueviolet", "#7B6E4F")]
pal <- colorRampPalette(brewer.pal(8,"Set1"))(8)
pal2 <- c("gray80", "dodgerblue", "goldenrod", "chartreuse", pal[c(1, 3:5, 8)], "blue", "gold1", "black")
length(pal2)

pdf(paste0("output/", date, "/", date, "_", cmpnd, "_JGI_genes_without_errorbars_normalized.pdf"), width = 16, height = 14)
pl <- ggplot(dat5, aes(x=time, y=mean, color=winners)) +
  geom_point() +
  labs(y = "nmol pNP", x = "Time (minutes)") +
  #geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
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

pdf(paste0("output/", date, "/", date, "_", cmpnd, "_JGI_genes_with_errorbars_normalized.pdf"), width = 16, height = 14)
pl <- ggplot(dat5, aes(x=time, y=mean, color=winners)) +
  geom_point() +
  labs(y = "nmol pNP", x = "Time (minutes)") +
  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
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

# Write to file
# write_csv(dat5, paste0("output/", date, "/", date, "_", cmpnd, "_JGI_normalized_data.csv"))

# Calculate slopes
a <- dat5 %>%
  ungroup() %>%
  dplyr::mutate(variable = as.character(variable)) %>%
  dplyr::filter(!grepl("Nothing|pNP", variable)) %>%
  dplyr::mutate(minutes = as.numeric(substr(time, 15, 16))) %>%
  dplyr::mutate(minutes = case_when(grepl(max(time), time) ~ 60,
                                    TRUE ~ minutes)) %>%
  dplyr::group_by(variable) %>%
  dplyr::slice(11:31) %>%
  dplyr::ungroup() %>%
  dplyr::select(variable, minutes, mean)


linr <- a %>%
  split(.$variable) %>%
  map(~lm(mean ~ minutes, data = .x)) %>%
  map_df(broom::tidy)

linr$variable <- sort(rep(unique(a$variable), 2))

wid <- dcast(linr, variable ~ term, value.var="estimate") %>%
  janitor::clean_names() %>%
  arrange(desc(minutes)) %>%
  dplyr::mutate(cmpnd = cmpnd)
wid
write_csv(wid, paste0("output/", date, "/", date, "_", cmpnd, "_slopes_calculated.csv"))

pdf(paste0("output/", date, "/", date, "_", cmpnd, "_slopes_plotted.pdf"), width = 20, height = 20)
    ggplot(a, aes(x = minutes, y = mean, color = variable)) +
      geom_point() +
      #  geom_smooth(method = "lm", se = TRUE)
      geom_abline(slope = wid$minutes, intercept = wid$intercept)
    # theme(legend.position="none")
dev.off()
