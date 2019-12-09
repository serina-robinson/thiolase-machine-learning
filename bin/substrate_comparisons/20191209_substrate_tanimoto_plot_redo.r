# Install packages
pacman::p_load("webchem", "readxl", "tidyverse", "ChemmineR", "RColorBrewer",
               "data.table", "plot3D", "plot3Drgl", "ggrepel")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the list of compounds
ll <- list.files("data/sigma_carboxylic_acid_scrape/", pattern = "C")
data_path <- "data/sigma_carboxylic_acid_scrape/"
rawdat <- tibble(filename = ll) %>%
  # purrr::map(read_excel) %>%   
  mutate(file_contents = map(filename,          # read files into
                             ~ read_excel(file.path(data_path, .))) # a new data column
  ) %>%
  unnest(.,cols = c(file_contents))  # read in all the files individually

dim(rawdat)
head(rawdat)

# Read in one plot
dat <- rawdat %>%
  janitor::clean_names() %>%
  dplyr::filter(!is.na(description)) %>%
  dplyr::mutate(cmpnd = case_when(grepl("acid", description) ~ paste0(word(description, sep = "acid", 1), "acid"),
                                  TRUE ~ trimws(paste0(word(description, sep = "[[:digit:]][[:digit:]]%", 1))))) %>%
  dplyr::mutate(cmpnd = gsub(" AldrichCPR", "", cmpnd)) %>%
  dplyr::mutate(cmpnd = gsub("\\(GC\\)", "", cmpnd)) %>%
  dplyr::mutate(cmpnd = gsub("\\(HPLC\\)", "", cmpnd)) %>%
  dplyr::mutate(cmpnd = gsub("≥[[:digit:]][[:digit:]]\\.[[:digit:]]%", "", cmpnd)) %>%
  dplyr::mutate(cmpnd = gsub("\\(~", "", cmpnd)) %>%
  dplyr::mutate(cmpnd = trimws(cmpnd)) %>%
  dplyr::distinct(cmpnd, .keep_all = TRUE) %>%
  dplyr::mutate(carbon_num = word(filename, 1, sep = "\\.xlsx"))
head(dat)

c_length <- dat %>%
  select(carbon_num, cmpnd)
head(c_length)
# dat$cmpnd[!grepl("acid", dat$cmpnd)]
# write_csv(data.frame(cbind(dat$description, dat$cmpnd)), "output/test_cmpnd_description.csv")


length(dat$cmpnd)
test_pull <- dat$cmpnd
table(duplicated(dat$cmpnd))

# all_cids <- get_cid(test_pull)
# nams <- names(all_cids)
# cids <- unlist(all_cids)
# 
# cid_df <- data.frame(as.matrix(all_cids)) %>%
#   rename(cids = as.matrix.all_cids.) %>%
#   unlist() %>%
#   na.omit() %>%
#   as.data.frame()
# colnames(cid_df)
# 
# cid_df$cmpnd <- rownames(cid_df)
# head(cid_df)
# colnames(cid_df)[1] <- "cid"
# final_df <- cid_df %>%
#   mutate(cid = as.numeric(cid))
# 
# key <- final_df
# head(key)
# key$cmpnd <- gsub("cids\\.", "", key$cmpnd)

#write_csv(key, "output/CID_cmpnd_name_key.csv")
# https://pubchem.ncbi.nlm.nih.gov/pc_fetch/pc_fetch.cgi
# write_delim(data.frame(final_df$cid), "data/CA_cids_output.txt", delim = "\t", col_names = F)

# Now pull CIDs for substrates that have been tested already or will be
sub_test <- read_excel("data/sigma_carboxylic_acid_scrape/20191209_final_substrates_tested.xlsx") %>%
  dplyr::pull(1)
sub_test
all_cids <- get_cid(sub_test)
nams <- names(all_cids)
nams
cids <- unlist(all_cids)
cids

cid_df <- data.frame(as.matrix(all_cids)) %>%
  rename(cids = as.matrix.all_cids.) %>%
  unlist() %>%
  na.omit() %>%
  as.data.frame()
colnames(cid_df)

cid_df$cmpnd <- rownames(cid_df)
head(cid_df)
colnames(cid_df)[1] <- "cid"
final_df <- cid_df %>%
  mutate(cid = as.numeric(cid))
final_df
# # 
key <- final_df
head(key)
key$cmpnd <- gsub("cids\\.", "", key$cmpnd)
key
# write_csv(key, "data/20191209_final_substrates_tested_CID_cmpnd_name_key.csv")

# write_delim(data.frame(final_df$cid), "data/20191209_final_substrates_tested_CA_cids_output.txt", delim = "\t", col_names = F)

# Combine everything
big_key <- read_csv("data/CID_cmpnd_name_key.csv")
small_key <- read_csv("data/20191209_final_substrates_tested_CID_cmpnd_name_key.csv")
merg_key <- big_key %>%
  bind_rows(small_key)
dim(merg_key)

# write_csv(merg_key, "data/20191209_final_1120_merged_CID_cmpnd_name_key.csv")
# cids_to_pull <- merg_key$cid
# write_delim(data.frame(cids_to_pull), "data/20191209_final_1120_merged_CID_cmpnd_name_key.txt", delim = "\t", col_names = F)
# 
# final_df <- read_delim("data/20191209_final_1120_merged_CID_cmpnd_name_key.txt", delim = "\t", col_names = F)

key <- read_csv("data/")
sdfset1 <- read.SDFset("data/20191209_final_substrates_tested_sdf")
sdfset2 <- read.SDFset("data/merged_sigma_CAs.sdf")

sdfset <- c(sdfset1, sdfset2)
# sdfset@ID <- makeUnique(sdfset@ID)
sdfset@ID[1:15] <- key$cmpnd

cid_pulled <- as.numeric(unlist(lapply(1:length(sdfset), function(x) { header(sdfset[[x]])[1] })))
length(unique(cid_pulled))


apset <- sdf2ap(sdfset)


# Identity
status <- read_excel("data/sigma_carboxylic_acid_scrape/20191209_final_substrates_tested.xlsx")
status  
key <- read_csv("data/20191209_final_substrates_tested_CID_cmpnd_name_key.csv")
dim(key)
# key$cmpnd[grep("chlorophenoxy", key$cmpnd)]

key_jnd <- key %>%
  left_join(., status, by = "cmpnd")
  # dplyr::mutate(status = case_when(is.na(status) ~ "commercially available",
  #                                  TRUE ~ status))

dups <- key[duplicated(key$cid),]
dups
newnams <- key[!duplicated(key$cid),]
newnams
newnamvec <- newnams


length(apset)
# clusters <- cmp.cluster(db=apset,  cutoff=0.7, save.distances = "data/distmat.rda")
# clusters <- cmp.cluster(db=apset,  cutoff = 1.0, save.distances = "data/distmat.rda")
head(clusters)
# write_csv(clusters, "output/sigma_CA_tanimoto_clusters.csv")
cluster.sizestat(clusters) # 56, 25 ,2


clusviz <- cluster.visualize(apset, clusters, size.cutoff = 1, quiet = TRUE)
x <- clusviz[,1]
y <- clusviz[,2]
trdat <- data.frame(cbind(x,y))


rownames(trdat) <- names(x)
tail(trdat)
# table(grepl("CMP", rownames(trdat)))
trdat$cmpnd <- rownames(trdat)

# head(trdat)
# dim(trdat)
# 
# tested_comps <- read_excel("data/sigma_carboxylic_acid_scrape/Larry_preferred_substrates.xlsx", col_names = F)


# Now try to find the substrates
subst <- sdf2ap(sdfset2)

which(newnams$cid %in% dups$cid[10])
for(i in 1:length(dups$cid)) {
  ind <- which(newnams$cid %in% dups$cid[i]) 
  newnamvec$cmpnd[ind] <- dups$cmpnd[i]
}
newnamvec$cmpnd
subst@ID <- newnamvec$cmpnd


subgrp <- trdat %>%
  left_join(., key_jnd, by = "cmpnd") %>%
  dplyr::mutate(cmpnd_labeled = case_when(is.na(status) ~ "NA",
                                          duplicated(cmpnd) ~ "NA",
                                          TRUE ~ cmpnd))
subgrp  

table(subgrp$cmpnd_labeled)
subgrp$cmpnd_labeled[!is.na(subgrp$status)]
subgrp$cmpnd_labeled <- gsub("NA", NA, subgrp$cmpnd_labeled)

#table(grepl(paste0(Larrys_comps$...1, collapse = "|"), subgrp$cmpnd))
# tofind <- subgrp$cmpnd[grepl(paste0(collapse = "|", Larrys_comps$...1[!Larrys_comps$...1 %in% subgrp$cmpnd]), subgrp$cmpnd)]


numseqs <- nrow(subgrp)
#pal2 <- c(colorRampPalette(brewer.pal(8, "Set2"))(4), "gray80")
pal2 <- c(colorRampPalette(brewer.pal(8, "Set1"))(8), "gray70")
pal2



pdf(paste0("output/", numseqs,"_PCoA_unlabeled_CA_Sigma.pdf"), width = 15, height = 15)
par(mar=c(0.01, 0.01, 0.01, 0.01))
ggplot(data = subgrp, aes(x = x, y = y, label = cmpnd_labeled)) + 
  scale_fill_manual(values = pal2) +
  scale_color_manual(values = pal2) +
  # geom_text(label = trdat$nms) +
  # geom_text(x = cdat$x, y = cdat$y, label = cdat$labl, check_overlap = T) +
  #geom_point(fill = "gray80") +
  geom_point(aes(fill = as.factor(status)), shape = 21, size = ifelse(is.na(subgrp$status), 4, 6)) +
  # scale_shape(solid = TRUE) +
  #geom_label_repel(force = 30, size = 4) +
  # labs(x=xlab,y=ylab) +
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.box = "horizontal",
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=14, face="bold", hjust=0))
dev.off()

pdf(paste0("output/", numseqs,"_PCoA_unlabeled_CA_Sigma.pdf"), width = 15, height = 15)
par(mar=c(0.01, 0.01, 0.01, 0.01))
ggplot(data = subgrp, aes(x = x, y = y, label = cmpnd_labeled)) + 
  scale_fill_manual(values = pal2) +
  scale_color_manual(values = pal2) +
  # geom_text(label = trdat$nms) +
  # geom_text(x = cdat$x, y = cdat$y, label = cdat$labl, check_overlap = T) +
  #geom_point(fill = "gray80") +
  geom_point(aes(fill = as.factor(status)), shape = 21, size = ifelse(is.na(subgrp$status), 4, 6)) +
  # scale_shape(solid = TRUE) +
 #  geom_label_repel(force = 30, size = 4) +
  # labs(x=xlab,y=ylab) +
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.box = "horizontal",
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=14, face="bold", hjust=0))
dev.off()
