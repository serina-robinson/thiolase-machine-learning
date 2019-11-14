# Pipeline for files with # Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer", 
               "ggplot2", "broom", "maditr")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

### Fill in your parameters of interest ###

### 2019-10-24
# cmpnd <- "Butoxy"
cmpnd <- "hexanoate"
date <- "2019-11-08"


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

resbind <- res[[1]] %>%
  bind_rows(res[[2]]) %>%
  dplyr::filter(!grepl("Nothing_", variable))

resbind <- resbind[complete.cases(resbind),] # Removes the "NOTHING's"

# Now normalize to the pNP standard curve

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

# Set random number seed for random color palette
set.seed(1234)
pal <- colorRampPalette(brewer.pal(8,"Set1"))(8)
pal2 <- c("gray80", "dodgerblue", "goldenrod", "chartreuse", pal[c(1, 3:5, 8)], "blue", "gold1", "black", distinctColorPalette(60))


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
  dplyr::slice(0:46) %>%
  dplyr::ungroup() %>%
  dplyr::select(variable, minutes, mean)


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
res[[1]]$org

names(res) <- orgs
resl <- plyr::ldply(res, data.frame)
resll <- do.call(rbind.data.frame, res)


# Find max slope for each organism
resmax <- resl %>%
  dplyr::filter(r2 >= 0.9) %>% # make sure R^2 is above or equal to 0.9
  group_by(org) %>%
  summarise_each(funs(max_slope = max), slope) %>%
  dplyr::filter(max_slope > 0)


# Merge with the original dataset
slope_merg <- resmax %>%
  inner_join(., resl, by = "org") %>%
  dplyr::filter(r2 >= 0.9) %>%
  group_by(org) %>%
  dplyr::filter(slope == max(slope)) %>%
  dplyr::select(org, max_slope, r2, intercept) %>%
  arrange(desc(max_slope))

# Plot winners on graph
merg_all <- slope_merg %>% 
  right_join(., a, by = c("org" = "variable")) %>%
  #left_join(., a, by = c("org" = "variable")) %>% # to exclude inactive ones
  dplyr::mutate(winners = case_when(is.na(max_slope) ~ " inactive",
                                    TRUE ~ org))

halo <- merg_all[grep("Halobacteriovorax", merg_all$org),]
halo
pdf(paste0("output/", date, "/", date, "_", cmpnd, "_JGI_genes_test_assessment_Halobacteriovorax_only.pdf"),  width = 16, height = 14)
pl <- ggplot(halo,  aes(x = minutes, y = mean, color = winners)) + 
  geom_point(alpha = ifelse(halo$winners == " inactive", 0.2, 1)) +
  geom_abline(slope = unique(halo$max_slope), intercept = unique(halo$intercept), color = pal2[1:length(unique(halo$max_slope))]) +
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
  scale_color_manual(values=pal2) 
# legend.position = "none")
pl
dev.off() 

# slope_sort <- slope_merg[order(slope_merg$max_slope, decreasing = T),]

pet_ind <- grep("Pet", slope_merg$org) - 1
pet_ind

# Visually assess results and remove any that look strange
# For dodecanoate, Microlunatus phosphovorus does not look active
cutoff <- 0
slope_final <- slope_merg[slope_merg$max_slope >= cutoff,]
if(length(pet_ind) > 0){
  slope_final <- slope_final[1:pet_ind,]
}
# slope_final <- slope_merg
# slope_final

write_csv(slope_merg, paste0("output/", date, "/", date, "_", cmpnd, "_all_data_calculated_slopes_round3.csv"))

# Plot updated winners on graph
merg_all_final <- slope_final %>% 
  right_join(., a, by = c("org" = "variable")) %>%
  #left_join(., a, by = c("org" = "variable")) %>% # to exclude inactive ones
  dplyr::mutate(winners = case_when(is.na(max_slope) ~ " inactive",
                                    TRUE ~ org))

pdf(paste0("output/", date, "/", date, "_", cmpnd, "_JGI_genes_linear_slopes_normalized_final_", cutoff, "_cutoff.pdf"), width = 20, height = 14)
pl <- ggplot(merg_all_final,  aes(x = minutes, y = mean, color = winners)) + 
  geom_point(alpha = ifelse(merg_all$winners == " inactive", 0.2, 1)) +
  geom_abline(slope = unique(merg_all_final$max_slope), intercept = unique(merg_all_final$intercept), color = pal2[1:length(unique(merg_all_final$max_slope))]) +
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
  scale_color_manual(values=pal2) 
# 
pl
dev.off()

