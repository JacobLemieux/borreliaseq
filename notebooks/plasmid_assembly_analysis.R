### Plasmid assembly analysis
# December 28 2022
# lemieux@broadinstitute.org

library(tidyverse)
library(Biostrings)

ref_genome <- readDNAStringSet("~/Dropbox/Tick/borrelia/genome/B_burgdorferi_B31.fa")
ref_lengths <- data.frame(Reference = sapply(names(ref_genome), function(x) strsplit(x, " ")[[1]][1]), Length = width(ref_genome))
metadata <- read_csv("~/Dropbox/Tick/borrelia/borreliaseq/results/LymeSeq_SampleTrack - Renaming Scheme Table for Jupyter.2023-01-03.csv") %>%
  filter(!is.na(Rename_A)) #%>%
  #filter(Location != "Slovenia") %>% 
  #filter(Rename_A != "UNY152")
replicon_db <- data.frame(Reference=NA)
for(i in 1:nrow(metadata)){
  contig_to_open <- metadata$Rename_A[i]
  contig <- read_tsv(paste("~/Dropbox/Tick/borrelia/borreliaseq/results/quast/contigs_reports/all_alignments_", contig_to_open, 
                           "_200.tsv", sep="")) %>%
    filter(Best_group == TRUE) %>% 
    mutate(E1 = as.numeric(E1), E2 = as.numeric(E2), S1 = as.numeric(S1), S2 = as.numeric(S2)) %>%
    mutate(ref_aligned = E1 - S1, contig_aligned = abs(E2 - S2))
  replicons <- contig %>% group_by(Reference) %>% 
    summarize(ref_align = sum(ref_aligned), contig_sum = sum(contig_aligned)) 
  names(replicons)[2:3] <- paste(names(replicons)[2:3], contig_to_open)
  replicons <- replicons %>% 
    select(-c(names(replicons[2])))
  replicon_db <- full_join(replicon_db, replicons, by = "Reference") %>% filter(!is.na(Reference))
}

replicon_db <- full_join(replicon_db, ref_lengths, by = "Reference")
names(replicon_db) <- gsub("contig_sum ", "", names(replicon_db))
replicon_PAM <- (replicon_db %>% select(-Reference)) / replicon_db$Length
replicon_binary <- ifelse(replicon_PAM > 0.5, 1, 0)
replicon_binary[is.na(replicon_binary)] <- 0
replicon_binary <- as_tibble(cbind(replicon_db$Reference, replicon_binary))
names(replicon_binary)[1] <- "Replicon"
replicon_binary <- data.frame(replicon_binary)
rownames(replicon_binary) <- replicon_binary$Replicon
replicon_binary <- replicon_binary %>% select(-Replicon)


pf32 <- data.frame(read_tsv("~/Dropbox/Tick/borrelia/borreliaseq/plasmid_types.tsv"))
rownames(pf32) <- pf32$Rename_A
pf32 <- pf32 %>% select(-Rename_A)
colnames(pf32) <- sapply(colnames(pf32), function(x) tail(strsplit(x, split = "_")[[1]],1))
colnames(pf32) <- gsub("ID.23.B31", "", colnames(pf32))
colnames(pf32) <- gsub("\\.", "-", colnames(pf32))
colnames(pf32)
pf32 <- data.frame(t(pf32))

colnames(pf32) <- paste("p", colnames(pf32), sep="")
pf32$Replicon <- rownames(pf32)
replicon_binary$Replicon <- rownames(replicon_binary)
plasmids <- left_join(replicon_binary, pf32, by = "Replicon")

write_tsv(replicon_binary, "~/Dropbox/Tick/borrelia/borreliaseq/plasmid_alignment_ref.tsv")
