### measure distance to PFam32 sequences
# Feb 6 2023
# lemieux@broadinstitute.org

library(Biostrings)
library(phangorn)
library(usedist)
library(tidyverse)

# functions

findPF32 <- function(testAAstring){
  D <- stringDist(c(testAAstring, casjens))
  X <- as.matrix(D)
  X <- X[-c(1),] / width(testAAstring)
  if(min(X[,1]) < 0.05){
    return(paste(rownames(X)[which(X[,1] == min(X[,1]))], sep="_"))
  }else{return("unknown")}
}

# load data
casjens <- readAAStringSet("~/Dropbox/Tick/borrelia/borreliaseq/blast/Casjens_pfam32_protein_reformatted.fa")

metadata <- read_csv("~/Dropbox/Tick/borrelia/borreliaseq/results/LymeSeq_SampleTrack - Renaming Scheme Table for Jupyter.2023-01-14.csv") %>%
  filter(!is.na(Rename_A)) 

# iterate through isolates, find closest-matching plasmid
plasmid_map <- data.frame(plasmid = names(casjens))
for(j in 1:length(metadata$Rename_A)){
#for(j in 1:5){
  PF32 <- readAAStringSet(paste("~/Dropbox/Tick/borrelia/borreliaseq/results/annotation/short_read2/pfam32/", metadata$Rename_A[j] ,"_pfam_subset.out", sep=""), format = "fasta")
  pf32_names <- data.frame(isolate = rep(1, length(PF32)), plasmid = NA)
  for(i in 1:length(PF32)){
   pf32_names$plasmid[i] <- findPF32(PF32[i])
   #print(findPF32(PF32[i]))
  }
  pf32_names <- pf32_names %>% distinct(isolate,plasmid)
  names(pf32_names)[1] <- metadata$Rename_A[j]
  plasmid_map <- full_join(plasmid_map, pf32_names, by = "plasmid")
}

# write out data
plasmid_map_out <- data.frame(plasmid_map)
rownames(plasmid_map_out) <- plasmid_map_out[,1]
plasmid_map_out <- plasmid_map_out[,-c(1)]
plasmid_map_out <- plasmid_map_out[!(rowSums(is.na(plasmid_map_out)) == ncol(plasmid_map_out)),]
plasmid_map_out <- t(plasmid_map_out)  
write_tsv(cbind(metadata$Rename_A, data.frame(plasmid_map_out)), "~/Dropbox/Tick/borrelia/borreliaseq/plasmid_types_v3.tsv")
