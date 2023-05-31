### Format genbank
# Jan 14 2023
# lemieux@broadinstitute.org

library(tidyverse)

genbank <- read_tsv("~/Dropbox/Tick/borrelia/borreliaseq/results/genbank/Pathogen.cl.1.0.tsv", skip = 11)

metadata <- read_csv("~/Dropbox/Tick/borrelia/borreliaseq/results/LymeSeq_SampleTrack - Renaming Scheme Table for Jupyter.2023-02-14.csv", na = c("", "NA", "ND")) %>% 
  filter(!is.na(Rename_A)) %>% 
  mutate(sample_name = Rename_A) %>% 
  mutate(sample_title = Original_Name) %>%
  mutate(organism = "Borrelia burgdorferi") %>% 
  mutate(collected_by = NA) %>% 
  mutate(collection_date = Year) %>%
  mutate(geo_loc_name = Origin) %>%
  mutate(geo_loc_name = recode(geo_loc_name, US = "USA", EU = "Slovenia")) %>%
  mutate(host = "Homo sapiens") %>% 
  mutate(host_disease = "Lyme disease") %>% 
  mutate(isolation_source = Source...12) %>% 
  mutate(lat_lon = NA) %>% 
  mutate(strain = Rename_A) 

metadata2 <- metadata %>% 
  select(c(sample_name, sample_title, organism, collected_by, collection_date, geo_loc_name, host, host_disease, isolation_source, lat_lon, strain ))

write_tsv(metadata2, "~/Dropbox/Tick/borrelia/borreliaseq/results/genbank/mod/Pathogen.cl.1.0_Lemieux_et_al_BBss.tsv")

# calculate coverage and output assembly metadata

x <- read_tsv("~/Dropbox/Tick/borrelia/borreliaseq/results/coverage_db.txt", col_names=FALSE) 
names(x) <- c("Name", "Number", "Length", "Coverage")

weighted_cov <- x %>% group_by(Name) %>% 
  summarize(weightedCov = sum(Coverage*Length)/sum(Length))

coverage <- x %>% group_by(Name) %>% 
  #filter(Length > 2000) %>%
  summarize(coverage = mean(Coverage))

y <- read_tsv("~/Dropbox/Tick/borrelia/borreliaseq/results/quast/report.tsv")
y <- data.frame(y)
rownames(y) <- y[,1]
y <- y[,-1]
y <- as.tibble(t(y))

# report numbers for the paper

summary(weighted_cov$weightedCov)
summary(coverage$coverage)


summary(as.numeric(y$N50))
summary(as.numeric(y$`Total length`))
summary(as.numeric(y$`# contigs`))


assem_stats <- inner_join(metadata, weighted_cov, by = c("sample_name" = "Name")) %>%
  mutate(genome_coverage = weightedCov) %>% 
  mutate(assembly_method = "SPAdes") %>% 
  mutate(assembly_method_version = "3.14.1") %>%
  mutate(sequencing_technologies = "Illumina") %>%
  mutate(filename = paste(sample_name, "_200.fa", sep="")) %>% 
  select(c(sample_name, assembly_method, assembly_method_version, sequencing_technologies, genome_coverage, filename))

write_tsv(assem_stats, "~/Dropbox/Tick/borrelia/borreliaseq/results/genbank/mod/genome_annotation.tsv")

# Generate ST 3

metadata <- inner_join(metadata, weighted_cov, by = c("Rename_A" = "Name")) %>% 
  mutate(Region = recode(Location, Nantucket = "Northeast (CT,RI,MI)", NY = "NY",
                         RI = "Northeast (CT,RI,MI)", CT = "Northeast (CT,RI,MI)", WI = "US Midwest",
                         Slovenia = "EU Slovenia"
  )) 
metadata3 <- metadata %>% 
  select(c(sample_name, sample_title, organism, collection_date, Region, host, host_disease, isolation_source, MLST, OspC_Type, RST_Type, `Multiple_EM(Y/N)`, `BC + (Y/N)`, `PCR + (Y/N)`, Disseminated )) %>%
  mutate(`PCR + (Y/N)` = recode(`PCR + (Y/N)`,"poz" = "Y", "neg" = "N")) %>%
  mutate(`BC + (Y/N)` = recode(`BC + (Y/N)`,"poz" = "Y", "neg" = "N"))
write_tsv(metadata3, "~/Dropbox/Tick/borrelia/borreliaseq/results/tables/submission/Supplemental_Table_2.tsv")

# figure out which proportion of isolates had either PCR or BCx info

bcx_done <- metadata3 %>% filter(`BC + (Y/N)` %in% c("Y" ,"N"))
pcr_done <- metadata3 %>% filter(`PCR + (Y/N)` %in% c("Y" ,"N"))
both_done <- union(bcx_done$sample_name, pcr_done$sample_name)

# assemble table summarizing metadata3 by location
metadata4 <- metadata3 %>% 
  mutate(BCx = recode(`BC + (Y/N)`, "Y" = "BCx pos", "N" = "BCx neg")) %>%
  mutate(PCR = recode(`PCR + (Y/N)`, "Y" = "PCR pos", "N" = "PCR neg")) %>%
  mutate(MEM = recode(`Multiple_EM(Y/N)`, "Y" = "MEM present", "N" = "MEM absent")) %>%
  mutate(From = recode(isolation_source, "Blood" = "Isolated from blood", "CSF" = "Isolated from CSF", "Skin" = "Isolated from skin"))
metadata_sum <- rbind(table(metadata4$Region),
      table(metadata4$BCx,metadata4$Region), 
      table(metadata4$PCR, metadata4$Region),
      table(metadata4$MEM,metadata4$Region),
      table(metadata4$From, metadata4$Region))
metadata_sum <- data.frame(metadata_sum, Totals = rowSums(metadata_sum))

single_EM <- metadata4 %>% filter(MEM == "MEM absent") %>%
  mutate(BCx = replace_na(BCx, "Single EM and BCx NA")) %>%
  mutate(PCR = replace_na(PCR, "Single EM and PCR NA")) %>% 
  mutate(BCx2 = recode(BCx, "BCx pos" = "Single EM and BCx pos", "BCx neg" = "Single EM and BCx neg")) %>%
  mutate(PCR2 = recode(PCR, "PCR pos" = "Single EM and PCR pos", "PCR neg" = "Single EM and PCR neg"))

metadata_sum2 <- rbind(table(single_EM$BCx2, single_EM$Region),
                       table(single_EM$PCR2, single_EM$Region))
metadata_sum2 <- data.frame(metadata_sum2, Totals = rowSums(metadata_sum2))


single_EM_sub <- single_EM %>% 
  filter(BCx == "Single EM and BCx NA" & PCR == "Single EM and PCR NA")

single_EM_sum <- t(data.frame(table(single_EM_sub$Region))) # there's got to be a better way to reformat this table...
colnames(single_EM_sum) <- single_EM_sum[1,]
single_EM_sum <- data.frame(single_EM_sum)
single_EM_sum <- single_EM_sum[-c(1),]
rownames(single_EM_sum) <- "Single_EM_with_PCR_and_BCx_NA"
single_EM_sum <- data.frame(single_EM_sum, Totals = sum(table(single_EM_sub$Region)))

metadata_sum_total <- rbind(metadata_sum, metadata_sum2, single_EM_sum)
rownames(metadata_sum_total)[1] <- "Total"
write.csv(metadata_sum_total, "~/Dropbox/Tick/borrelia/borreliaseq/results/tables/Submission/Supplemental_Table_1.csv")      
