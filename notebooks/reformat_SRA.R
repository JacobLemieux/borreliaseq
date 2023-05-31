### SRA submission
# lemieux@broadinstitute.org
# May 25 2023

library(xmlconvert)
library(tidyverse)
xml <- xml_to_df("../results/genbank/biosample_result.xml", records.tags = "BioSample") 
xml2 <- xml %>% 
  select(Ids) %>%
  separate(Ids, sep="\\|\\|", into = c("BioSample", "ID")) %>%
  mutate(BioSample = gsub("\\|", "", BioSample)) %>%
  mutate(BioSample = gsub("Id~", "", BioSample)) %>%
  mutate(ID = gsub("\\|", "", ID)) %>%
  mutate(ID = gsub("Id~", "", ID)) %>% 
  mutate(R1 = paste(ID, "R1.fastq.gz", sep="_")) %>%
  mutate(R2 = paste(ID, "R2.fastq.gz", sep="_"))

write_csv(xml2,"../results/genbank/sample_accession_ID.csv")
