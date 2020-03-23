extract_12angstrom_one_hot <- function(query_fil) {

  library(seqinr)
  library(bgafun)
  
  # Read in the reference sequence
  query <- readAAStringSet(query_fil)
  ref <- readAAStringSet("data/alignments/4KU5.fasta")
  names(ref) <- "4KU5"
  
  # Align the query and the reference
  alned <- AlignSeqs(c(ref, query), verbose = TRUE)
  # alned <- AAStringSet(muscle(c(ref, query), in1 = "data/A_domains_muscle.fasta", in2 = query_fil, profile = T)) # If you want to align with a profile
  query_aln <- alned[length(alned)]
  ref_aln <- alned["4KU5"]
  
  # Read in the 12 angstrom residues
  aa84_inds <- read_csv("data/machine_learning/84_residues_12_angstroms_4KU5_S143.csv") %>%
    pull()
  aa84_inds_adj <- aa84_inds # Depends on ref

  # Exract the 34 amino acid positions
  poslist <- list()
  position = 1
  
  for(i in 1:width(ref_aln)) {
    if (substr(ref_aln, i, i) != "-") {
      if (position %in% aa84_inds_adj) {
        poslist[[i]] <- i
      }
    position = position + 1
    }
  }
  
  # Get the new indices
  new_84inds <- unlist(poslist)
  new_84inds
  
  # Get 84 aa code
  query_pos <- as.character(unlist(lapply(1:length(new_84inds), function(x) {
    substr(query_aln, new_84inds[x], new_84inds[x]) })))
  return(query_pos)
}



  