extract_12angstrom_aas <- function(query_fil) {

  library(seqinr)
  library(bgafun)
  
  # Read in the reference sequence
  query <- readAAStringSet(query_fil)
  ref <- readAAStringSet("data/4KU5.fasta")
  names(ref) <- "4KU5"
  
  # Align the query and the reference
  alned <- AlignSeqs(c(ref, query), verbose = FALSE)
  query_aln <- alned[length(alned)]
  ref_aln <- alned["4KU5"]
  
  # Read in the 12 angstrom residues
  aa84_inds <- read_csv("data/84_residues_12_angstroms_4KU5_S143.csv") %>%
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
  query_84 <- AAStringSet(paste0(query_pos, collapse = ""))
  writeXStringSet(query_84, "output/temp_84_aa.fasta")
  
  # Convert the 713 aa signatures to features
  rdaln <- read.alignment(file = 'output/temp_84_aa.fasta', format = "fasta")
  rdaln$seq <- toupper(rdaln$seq)
  aa <- bgafun::convert_aln_AAP(rdaln) #5 physicochemical properties
  
  aadf <- data.frame(aa, stringsAsFactors = F)
  head(aadf)
  
  aap <- aadf %>%
    dplyr::mutate(nms = rownames(.)) %>%
    dplyr::select(-contains("D")) %>%
    dplyr::filter(!grepl("ERROR", X1A)) %>%
    dplyr::select(-nms)
  colnames(aap) <- gsub("^X","",colnames(aap))
  colnames(aap) <- paste0(c("polrty", "secstr", "molsz", "elechrg"), "_", colnames(aap), "_4KU5pos_", rep(aa84_inds, each = 4))
  
  numfeats <- length(colnames(aap)) # 1096

  feats <- aap

  return(feats)
}



  