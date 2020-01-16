## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'data.table', 'tidyverse')

## Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the pid matrix
pid <- fread("data/alignments/73_JGI_OleA_wildtype_Xc_pid.txt", skip = 1, sep = " ", data.table = F)
pid

pid_num <- pid %>%
  select(-V1, -V2)
rowMeans(pid_num)
summary(na.omit(apply(pid_num, 1, min)))
which.min(na.omit(apply(pid_num, 1, min)))

na.omit(apply(pid_num, 1, min))

tmp <- na.omit(apply(pid_num, 1, min))
n <- length(tmp)
tmp
sort(tmp, partial = n - 1)[n-1]

min( tmp[tmp!=0.00] ) # 2.78% 
