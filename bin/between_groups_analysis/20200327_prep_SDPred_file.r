# Install packages
pacman::p_load(DECIPHER, tidyverse, seqinr, bgafun, RColorBrewer)

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the dataset
dat <- readAAStringSet("data/machine_learning/40_OleA_broad_narrow_sub_spec_12angstrom.faa")
head(dat)
