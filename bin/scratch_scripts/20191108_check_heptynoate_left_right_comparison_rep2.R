# Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer", 
               "ggplot2", "data.table", "maditr", "broom")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

### 2019-11-08
date <- "2019-11-08"
cmpnd <- "heptynoate"

# Read in the heptynoate slopes
hept <- read_csv('output/2019-11-08/2019-11-08_heptynoate_biorep2_all_data_calculated_slopes.csv') %>%
  dplyr::mutate(genus = word(org, 1 , sep = " "))
table(hept$genus)

# Merge with the key
key <- read_csv('data/72_OleA_masterwell_org_key.csv')
key$orgs <- gsub("_", " ", key$orgs)
key$orgs <- gsub("Lysobacter", "Luteimonas", key$orgs)
keygen <- key %>%
  dplyr::mutate(genus = word(orgs, 1, sep = " "))
head(key)

keyjn <- hept %>%
  left_join(., keygen, by ="genus")
keyjn$master_well[keyjn$genus == "XC"] <- "XC"
# hept[grep("2-E9|1-A4|1-D2|1-E3|1-D2|1-G3|1-F3|1-F4|1-F5|1-B6|1-H5|1-A6", ),]
hept[grep("Luteimonas", hept$genus),]
key[grep("2-E9|1-A4|1-D2|1-E3|1-D2|1-G3|1-F3|1-F4|1-F5|1-B6|1-H5|1-A6", key$master_well),]

# Read in the plate template
temp <- read_excel(paste0("data/", date, "/", "Plate Set up SynBio Paper.xlsx"), skip = 3) %>%
  janitor::clean_names() %>%
  dplyr::select(-x1, -x12) %>%
  as.data.frame()#%>% #%>% # NOTE THAT X11 WILL BE FILLED 
#  as.matrix()
  # #t %>%
  # as.vector()
temp
tail(keyjn)

# Substitue master wells for activities
temp_org <- temp
for(i in 1:ncol(temp_org)) {
  for(j in 1:length(keyjn$master_well)) {
    ind <- grep(keyjn$master_well[j], temp[,i])
    if(length(ind) > 0){
      if(!is.na(ind)) {
        temp_org[ind, i] <- round(keyjn$max_slope[j], 2)
      }
    }
  }
}

write_csv(temp_org, "output/scratch_output/left_right_bio_rep_2_heptynoate_comparison.csv")

