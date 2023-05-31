### Main figures, revised manuscript
# April 28 2023
# lemieux@broadinstitute.org

# requires local path as borreliaseq/notebooks

# Read in libraries ####
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(tidyverse)
library(ape)
library(phytools)
library(cowplot)
source("../scripts/Heatmap_PAM.R")
library(GGally)
library("ggpubr")
#library(DECIPHER)
library(forcats)
library("DescTools")
library(ggpubr)
library(ggfortify)
library(ggnewscale)
library(ggstar)
library(scales)
library(ggrepel)
#library(beastio)

# Read in data ####
# read in trees and midpoint root
x <- read.newick("../results/roary-results/core_gene_alignment.nwk")
x1 <- read.newick("../results/cfOUT_.labelled_tree.newick")
y <- read.newick("../results/OspC_Bbss_combined_4_24_2021.nwk")
z <- read.newick("../results/roary-results/split/core_gene_alignment_split.tree")
a <- read.newick("../results/roary-results/accessory_binary_genes.fa.newick")
x <- midpoint.root(x)
x1 <- midpoint.root(x1)
y <- midpoint.root(y)
z <- midpoint.root(z)
a <- midpoint.root(a)
b <- read.beast("../results/BEAST/no_dates_ancestral_geography/core_gene_alignment.aln.mcc.tree")
iq_tree <- read.iqtree("../results/roary-results/core_gene_alignment.aln.contree")
iq_tree@phylo <- midpoint.root(iq_tree@phylo)
O_tree <- read.beast("../results/trees/OspC_MCC.tree")
# check and parse labels

# read in OspC tree and label tips
OspC_names <- O_tree@phylo$tip.label
OspC_names_table <- data.frame(str_split(OspC_names, pattern="_", simplify=TRUE))
dup_vector <- OspC_names_table[which(duplicated(OspC_names_table[,1])),1]
OspC_names_table$duplicated <- OspC_names_table[,1] == dup_vector
OspC_names_table$new_name <- ifelse(!OspC_names_table$duplicated, OspC_names_table[,1], paste(OspC_names_table[,1], OspC_names_table[,2], sep="_"))
O_tree@phylo$tip.label <- OspC_names_table$new_name

# read in metadata
metadata <- read_csv("../results/LymeSeq_SampleTrack - Renaming Scheme Table for Jupyter.2023-02-14.csv", na = c("", "NA", "ND")) %>% 
                   mutate(RST_Type = factor(RST_Type)) %>%
                   mutate(OspCK = ifelse(OspC_Type == "K", "K", "other")) %>% 
                   mutate(OspCA = ifelse(OspC_Type == "A", 1, 0)) %>% 
                   mutate(Region = recode(Location, Nantucket = "US Northeast", NY = "US Northeast",
                                          RI = "US Northeast", CT = "US Northeast", WI = "US Midwest",
                                          Slovenia = "EU Slovenia"
                   )) %>%
                   mutate(Severity_Score = recode(Severity, mild = 0, moderate = 1, severe = 2))
names(metadata)[1] <- "label"
metadata <- metadata %>%
  filter(label %in% x$tip.label)

# read in distances and infer clusters
d <- data.frame(read_tsv("../results/kWIP-results/Bb.dist"))
# clean up
rownames(d) <- d[,1]
d <- d[,-1]
d <- d[rownames(d) %in% metadata$label, colnames(d) %in% metadata$label]
#plot(cmdscale(d))
set.seed(1111)
clust <- kmeans(cmdscale(d), 3)$cluster %>%
  as.factor()
set.seed(1111)
clust2 <- kmeans(cmdscale(d), 4)$cluster %>%
  as.factor()
mds <- cmdscale(d) %>%
  as_tibble %>%
  mutate(groups = clust, label = rownames(cmdscale(d))) %>% 
  mutate(groups2 = clust2) %>%
  mutate(wgs_group = factor(recode(groups2, `1` = "A", `2`="B.2", `3`="B.1", `4`="C"))) %>%
  mutate(wgs_group = factor(wgs_group, levels(wgs_group)[c(1,3,2,4)])) %>% 
  mutate(WGS_Group = factor(recode(groups, `1` = "A", `2`="C", `3`="B"))) %>% 
  mutate(WGS_Group = factor(WGS_Group, levels(WGS_Group)[c(1,3,2)]))
# merge in MDS assignments to metadata
metadata_extended <- full_join(metadata, mds, by = "label")  %>% 
  mutate(dissem_bin = ifelse(Disseminated == "D", 1, 0)) %>% 
  mutate(Country = recode(Region, `US Northeast` = "US", `US Midwest` = "US", `EU Slovenia` = "Slovenia"))
metadata_extended <- data.frame(metadata_extended)

# Create basic trees ####
p <- ggtree(x)
p1 <- ggtree(y)

# Calculate proportions of RST and OspC types
print("RST Types")
table(metadata$RST_Type)/sum(table(metadata$RST_Type))
# Ensure that we have a full set of isolates
print("total number of isolates with RST typing")
sum(table(metadata[,5]))
print("OspC Types")
table(metadata$OspC_Type)/sum(table(metadata$OspC_Type))
print("total number of isolates with OspC typing")
sum(table(metadata$OspC_Type))

# Count numbers of OspC Types
table(metadata$OspC_Type)
sum(table(metadata$OspC_Type))

# Make a variety of trees
location_tree <- ggtree(x) %<+% metadata + 
  geom_tippoint(aes(color=Region, subset=!is.na(Location)), size = 4)

loc_ospC <- ggtree(y) %<+% metadata + 
  geom_tippoint(aes(color=Region, subset=!is.na(Location)), size = 4)

OspC_type_tree <- ggtree(x) %<+% metadata + 
  geom_tippoint(aes(color=OspC_Type, subset=!is.na(OspC_Type)), size = 4)

ospC_tree <- ggtree(y) %<+% metadata + 
  geom_tippoint(aes(color=OspC_Type, subset=!is.na(OspC_Type)), size = 4)

RST_type_tree <- ggtree(x) %<+% metadata + 
  geom_tippoint(aes(color=RST_Type, subset=!is.na(RST_Type)), size = 4)

RST_ospC <- ggtree(y) %<+% metadata + 
  geom_tippoint(aes(color=RST_Type, subset=!is.na(RST_Type)), size = 4)

dissem_tree <- ggtree(x) %<+% metadata_extended + 
  geom_tippoint(aes(color=Disseminated, subset=!is.na(Disseminated)), size = 4, alpha = 0.6)

source_tree <- ggtree(x) %<+% metadata + 
  geom_tippoint(aes(color=Source, subset=!is.na(Source)), size = 4)

MEM_tree <- ggtree(x) %<+% metadata + 
  geom_tippoint(aes(color=`Multiple_EM(Y/N)`, subset=!is.na(`Multiple_EM(Y/N)`)), size = 4)

MLST_tree <- ggtree(x) %<+% metadata + 
  geom_tippoint(aes(color=MLST, subset=!is.na(MLST)), size = 4)

MLST_ospC <- ggtree(y) %<+% metadata + 
  geom_tippoint(aes(color=MLST, subset=!is.na(MLST)), size = 4)

wgs_tree <- ggtree(x) %<+% metadata_extended + 
  geom_tippoint(aes(color=WGS_Group, subset=!is.na(groups)), size = 4)

# Annotate trees with WGS groups
wgs_tree_split <- ggtree(z) %<+% metadata_extended + 
  geom_tippoint(aes(color=groups, subset=!is.na(groups)), size = 4)

wgs_tree_RST <- ggtree(x) %<+% metadata_extended + 
  geom_tippoint(aes(color=RST_Type), size = 2)

ospc_tree <- ggtree(y) %<+% metadata_extended + 
  geom_tippoint(aes(color=groups, subset=!is.na(groups)), size = 4)

wgs2_tree <- ggtree(x) %<+% metadata_extended + 
  geom_tippoint(aes(color=wgs_group, subset=!is.na(wgs_group)), size = 4)
wgs_tree_ospc <- ggtree(x) %<+% metadata + geom_tippoint(aes(color = OspC_Type), size = 4)

# annotate tree with roary data and annotations
genotype <- as.data.frame(read_csv("../results/roary-results/gene_presence_absence.csv"))
rownames(genotype) <- genotype[,1]
genotype <- genotype[,-c(1:14)]
# convert to Rtab object
gene_mat <- mapply(function(x) ifelse(!is.na(x), 1, 0), genotype)
gene_mat <- data.frame(gene_mat)
gene_mat$Gene <- rownames(genotype)
# read in annotation from roary
annotation <- read.csv("../results/roary-results/gene_presence_absence.csv")
annotation <- annotation[,c(1,3)]
annotation$Gene_trunc <- gsub("group", "g", annotation$Gene)
annotation$Annotation <- gsub("hypothetical protein", "hyp", annotation$Annotation)

gt.annotated <- merge(annotation, gene_mat)
rownames(gt.annotated) <- paste(gt.annotated$Gene_trunc,":  ", gt.annotated$Annotation, sep="")
gt.annotated <- gt.annotated[,-c(1:3)]

# read in annotation from the Make_Manhattan.ipynb
annot <- as_tibble(read_csv("../results/annotation/compiled_annotation.csv")) %>% 
  mutate(Gene = group) %>% select(-group) %>% select(Gene,ID, Locus, B31_Annotation, Localization) %>% distinct()
gt.annot <- merge(annot, gene_mat)  #%>% filter(Localization == "S")
gt.annot.L <- merge(annot, gene_mat) %>% filter(Localization == "S")
rownames(gt.annot) <- paste(gt.annot$Locus, gt.annot$B31_Annotation, gt.annot$Gene)
rownames(gt.annot.L) <- paste(gt.annot.L$Locus, gt.annot.L$B31_Annotation, gt.annot.L$Gene)
gt.annot <- gt.annot[,-c(1:5)]
gt.annot.L <- gt.annot.L[,-c(1:5)]
gt.annot.E <- gt.annot.L[grep("erp",rownames(gt.annot.L), ignore.case=TRUE),]
gt.annot.M <- gt.annot.L[grep("mlp",rownames(gt.annot.L), ignore.case=TRUE),]

# compute number of ORFS; might be nice to be able to plot these with a barplot next to the tree, but haven't figured out how to do this yet.
genome_length <- rowSums(t(gt.annotated))
genome_lengths <- data.frame(ORF_length = genome_length, label = colnames(gt.annotated))
metadata_extended <- left_join(metadata_extended, genome_lengths, by = "label")


## 
x2 <- read.tree("../results/RAxML_bestTree.raxOUT.nwk")
x2 <- midpoint.root(x2)

## 
mds_plot2 <- ggscatter(mds, x = "V1", y = "V2", 
                       #label = rownames(d),
                       color = "wgs_group",
                       size = 1, 
                       ellipse = TRUE,
                       #ellipse.type = "convex",
                       repel = TRUE) + labs(x = "MDS Coordinate 1", y = "MDS Coordinate 2")
mds_plot2


# Generate MDS plots annotated with metadata ####

mds_OspC <- ggscatter(metadata_extended %>% filter(!is.na(OspC_Type)), x = "V1", y = "V2", 
                      #label = rownames(d),
                      color = "OspC_Type",
                      size = 1, 
                      repel = TRUE) + labs(x = "MDS Coordinate 1", y = "MDS Coordinate 2")

mds_OspC_lab <- ggscatter(metadata_extended %>% filter(!is.na(OspC_Type)), x = "V1", y = "V2", 
                      #label = rownames(d),
                      color = "OspC_Type", label = "OspC_Type",
                      size = 1, 
                      repel = TRUE) + labs(x = "MDS Coordinate 1", y = "MDS Coordinate 2")

mds_MLST <- ggscatter(metadata_extended %>% filter(!is.na(OspC_Type)), x = "V1", y = "V2", 
                      #label = rownames(d),
                      color = "MLST",
                      size = 1, 
                      repel = TRUE) + labs(x = "MDS Coordinate 1", y = "MDS Coordinate 2", title = "MLST") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

dissem_data <- ggscatter(metadata_extended, x = "V1", y = "V2", 
                         #label = rownames(d),
                         color = "Disseminated",
                         size = 1, 
                         repel = TRUE) + labs(x = "MDS Coordinate 1", y = "MDS Coordinate 2")

loc_data <- ggscatter(metadata_extended, x = "V1", y = "V2", 
                      #label = rownames(d),
                      color = "Region",
                      size = 1, 
                      repel = TRUE) + labs(x = "MDS Coordinate 1", y = "MDS Coordinate 2")
plot_grid(mds_OspC, mds_MLST, dissem_data, loc_data, labels = c("A", "B", "C", "D"))
#ggsave("../results/figures/paper/revision/Figure_S1.jpg", height = 10, width = 10)

mds_OspC


# Create core + accessory genome PAM ####

wgs_tree_ospc <- wgs_tree_ospc + 
  geom_treescale(x = 0.002, y = -25, width = 0.002, fontsize = 6, offset = 2)

upper_threshold = 0 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 10 # remove from plotting anything in less than lower_threshold of isolates
genotype.trunc <- gt.annotated[rowSums(gt.annotated)> lower_threshold & (rowSums(gt.annotated) < (ncol(gt.annotated)-upper_threshold)),]
p_hm <- heatMapPAM(wgs_tree,t(genotype.trunc), colnames=FALSE, col_colours="blue", 
                   colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

wgs_tree_reg <- ggtree(x) %<+% metadata + geom_tippoint(aes(color = Region))
p_hm_reg <- heatMapPAM(wgs_tree_reg,t(genotype.trunc), colnames=FALSE, col_colours="blue", colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)

wgs_tree_ospc <- ggtree(x) %<+% metadata + geom_tippoint(aes(color = OspC_Type), size = 2) 
  
p_hm_ospc <- heatMapPAM(wgs_tree_ospc,t(genotype.trunc), colnames=FALSE, col_colours="blue", colnames_angle = -90,
                        hjust =1, width=7, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)
ggsave("../results/figures/paper/revision/Figure_5A.jpg", height = 4.5, width = 8)


# Create core lipoproteome PAM #### 
upper_threshold = 20 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 10 # remove from plotting anything in less than lower_threshold of isolates
genotype.trunc2 <- gt.annot[rowSums(gt.annot)> lower_threshold & (rowSums(gt.annot) < (ncol(gt.annot)-upper_threshold)),]

wgs_tree_ospc <- ggtree(x) %<+% metadata + geom_tippoint(aes(color = OspC_Type), size = 2) 
  

p_hm2 <- heatMapPAM(wgs_tree_ospc,t(genotype.trunc2), colnames=TRUE, col_colours="blue", 
                    colnames_angle = 90,hjust =0, width=35, font.size=5, cluster_cols=TRUE, 
                    null_colour="white",colnames_offset_y=-20, border_colour=NULL,
                    colnames_position = "bottom")
#ggsave("../results/figures/core_genome_clusters_annotation2.pdf", height = 100, width = 200, limitsize=FALSE)
p_hm2 <- heatMapPAM(wgs_tree_reg,t(genotype.trunc2), colnames=TRUE, col_colours="blue", 
                    colnames_angle = 90,hjust =0, width=35, font.size=5, cluster_cols=TRUE, 
                    null_colour="white",colnames_offset_y=-20, border_colour=NULL,
                    colnames_position = "bottom")
ggsave("../results/figures/paper/revision/annotated_cluster_map.pdf", height = 100, width = 200, limitsize=FALSE)


# Create lipoproteome PAM #### 

# core lipoproteome
upper_threshold = 0 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 220 # remove from plotting anything in less than lower_threshold of isolates
genotype.L.trunc2 <- gt.annot.L[rowSums(gt.annot.L)>= lower_threshold & (rowSums(gt.annot.L) <= (ncol(gt.annot.L)-upper_threshold)),]

wgs_tree_ospc_lipo <- ggtree(x) %<+% metadata + geom_tippoint(aes(color = OspC_Type), size = 5) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "left", 
                                             title.theme = element_text(angle = 90, size = 20),
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 90, size = 20), override.aes = list(size=10)) ) + 
  geom_treescale(x = 0.002, y = -25, offset = 2, fontsize = 16, width = 0.002, linesize = 2)
wgs_tree_reg_lipo <- ggtree(x) %<+% metadata + geom_tippoint(aes(color = Region), size = 5) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "left", title.theme = element_text(angle = 90),
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 90, size = 20)))

p_hm2.L <- heatMapPAMlipo(wgs_tree_ospc_lipo,t(genotype.L.trunc2), colnames=TRUE, col_colours="blue", 
                      colnames_angle = 90,hjust =0, width=10, font.size=5, cluster_cols=TRUE, 
                      null_colour="white",colnames_offset_y=-180, border_colour=NULL,
                      colnames_position = "bottom")
ggsave("../results/figures/paper/revision/annotated_core_lipoproteom_cluster_map.pdf", height = 20, width = 13.5, limitsize=FALSE)
ggsave("../results/figures/paper/revision/Figure_6A.jpg", height = 20, width = 13.5, limitsize=FALSE)


# strain-variable/accessory

wgs_tree_ospc_lipo2 <- ggtree(x) %<+% metadata + geom_tippoint(aes(color = OspC_Type), size = 2) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "left", title.theme = element_text(angle = 90),
                                               label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                               label.theme = element_text(angle = 90, size = 20), override.aes = list(size=10))) + 
  geom_treescale(x = 0.002, y = -25, offset = 2, fontsize = 8, width = 0.002, linesize = 1.5)

upper_threshold = 56 # remove from plotting anything in more than N - upper_threshold of the isolates, i.e. plot only ORFs in 90+ % of isolates
lower_threshold = 15 # remove from plotting anything in less than lower_threshold of isolates
genotype.L.trunc2 <- gt.annot.L[rowSums(gt.annot.L) >= lower_threshold & (rowSums(gt.annot.L) < (ncol(gt.annot.L)-upper_threshold)),]

wgs_tree_ospc_lipo <- ggtree(x) %<+% metadata + geom_tippoint(aes(color = OspC_Type), size = 18) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "left", title.theme = element_text(angle = 90),
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 90)))
p_hm2.L <- heatMapPAMlipo(wgs_tree_ospc_lipo2,t(genotype.L.trunc2), colnames=TRUE, col_colours="blue", 
                    colnames_angle = 90,hjust =0, width=10, font.size=3, cluster_cols=TRUE, 
                    null_colour="white",colnames_offset_y=-180, border_colour=NULL,
                    colnames_position = "bottom", show_legend = FALSE)
ggsave("../results/figures/paper/revision/Figure_6B.jpg", height = 10, width = 16.5, limitsize=FALSE)
ggsave("../results/figures/paper/revision/annotated_accessory_lipoproteom_cluster_map.pdf", height = 10, width = 16.5, limitsize=FALSE)

# Erp / Mlp map
wgs2_tree <- ggtree(x) %<+% metadata_extended + 
  geom_tippoint(aes(color=wgs_group, subset=!is.na(wgs_group)), size = 4) + 
  scale_colour_discrete(name = "WGS Group",guide = guide_legend(direction = "horizontal", title.position = "left", 
                                            title.theme = element_text(angle = 90, size = 20),
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 90, size = 20))) + 
  geom_treescale(x = 0.002, y = -15, offset = 2, fontsize = 4, width = 0.002, linesize = 1)

p_hm2.E <- heatMapPAMlipo(wgs2_tree,t(gt.annot.E), colnames=TRUE, col_colours="blue", 
                      colnames_angle = 90,hjust =0, width=10, font.size=3, cluster_cols=TRUE, 
                      null_colour="white",colnames_offset_y=-100, border_colour=NULL,
                      colnames_position = "bottom")

ggsave("../results/figures/paper/revision/Figure_6C.jpg", height = 8, width = 8, limitsize=FALSE)

p_hm2.M <- heatMapPAMlipo(wgs2_tree,t(gt.annot.M), colnames=TRUE, col_colours="blue", 
                      colnames_angle = 90,hjust =0, width=10, font.size=3, cluster_cols=TRUE, 
                      null_colour="white",colnames_offset_y=-140, border_colour=NULL,
                      colnames_position = "bottom") 

ggsave("../results/figures/paper/revision/Figure_6D.jpg", height = 8, width = 10, limitsize=FALSE)

# calculate some genome stats as comparator to lipo stats

gt.annot.2 <- data.frame(label = colnames(gt.annot), orf_num = colSums(gt.annot))
genome_stats <- merge(gt.annot.2, metadata_extended %>% select(label, wgs_group, Region, dissem_bin, OspC_Type))
p6.0c <- genome_stats %>% ggplot(aes(x = orf_num, y = dissem_bin)) + geom_jitter(height = 0.01, alpha = 0.3) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial")) + 
  theme_bw() + 
  labs(x = "Number of ORFs", y = "Probability") + 
  coord_cartesian(x = c(1000, 1250))

genome_glm <- glm(dissem_bin ~ orf_num, family = "binomial", data = genome_stats)
summary(genome_glm)


# calculate the number of Erps, Mlps by OspC type
gt.annot.L2 <- data.frame(label = colnames(gt.annot.L), lipo_num = colSums(gt.annot.L))
lipo_stats <- merge(gt.annot.L2, metadata_extended %>% select(label, wgs_group, Region, dissem_bin, OspC_Type))
lipo_stats$OspC_Type <- as.factor(lipo_stats$OspC_Type)
p6.1c <- lipo_stats %>% ggplot(aes(x = lipo_num, y = dissem_bin)) + geom_jitter(height = 0.01, alpha = 0.3) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial")) + 
  theme_bw() + 
  labs(x = "Number of lipoproteins", y = "Probability")

lipo_glm <- glm(dissem_bin ~ lipo_num, family = "binomial", data = lipo_stats)
summary(lipo_glm)

lipo_comp = list(c("A", "B.1"), c("A", "B.2"), c("A", "C"), c("B.1", "B.2"), c("B.1", "C"), c("B.2", "C"))
p6.1 <- lipo_stats %>% ggstripchart(x = "wgs_group", y = "lipo_num") + 
  #facet_grid(~Region) +
  stat_compare_means(comparisons = lipo_comp, label = "p.signif") + 
  labs(x = "WGS Group", y = "Number of surface-exposed lipoproteins")

gt.annot.E2 <- data.frame(label = colnames(gt.annot.E), erp_num = colSums(gt.annot.E))
erp_stats <- merge(gt.annot.E2, metadata_extended %>% select(label, wgs_group, Region, OspC_Type, dissem_bin))
erp_stats$OspC_Type <- as.factor(erp_stats$OspC_Type)
p6.2c <- erp_stats %>% ggplot(aes(x = erp_num, y = dissem_bin)) + geom_jitter(height = 0.01, alpha = 0.3) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial")) +
  theme_bw() + 
  labs(x = "Number of Erps", y = "Probability")


erp_glm <- glm(dissem_bin ~ erp_num, family = "binomial", data = erp_stats)
summary(erp_glm)


p6.2 <- erp_stats %>% ggstripchart(x = "wgs_group", y = "erp_num") + 
  #facet_grid(~Region) +
  stat_compare_means(comparisons = lipo_comp, label = "p.signif") + 
  labs(x = "WGS Group", y = "Number of Erps")

gt.annot.M2 <- data.frame(label = colnames(gt.annot.M), mlp_num = colSums(gt.annot.M))
mlp_stats <- merge(gt.annot.M2, metadata_extended %>% select(label, wgs_group, Region, OspC_Type, dissem_bin))
mlp_stats$OspC_Type <- as.factor(mlp_stats$OspC_Type)

p6.3c <- mlp_stats %>% ggplot(aes(x = mlp_num, y = dissem_bin)) + geom_jitter(height = 0.01, alpha = 0.3) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial")) +
  theme_bw() + 
  labs(x = "Number of Mlps", y = "Probability")

mlp_glm <- glm(dissem_bin ~ mlp_num, family = "binomial", data = mlp_stats)
summary(mlp_glm)

p6.3 <- mlp_stats %>% ggstripchart(x = "wgs_group", y = "mlp_num") + 
  #facet_grid(~Region) +
  stat_compare_means(comparisons = lipo_comp, label = "p.signif") + 
  labs(x = "WGS Group", y = "Number of Mlps")

plot_grid(p6.1, p6.2, p6.3, nrow = 1)
ggsave("../results/figures/paper/revision/Figure_6C.jpg", height = 8, width = 8)

p6.1b <- lipo_stats %>% ggstripchart(x = "OspC_Type", y = "lipo_num", color = "wgs_group") + 
  #facet_grid(~Region) +
  labs(x = "OspC Type", y = "Lipoproteins") + 
  guides(color = guide_legend(title = "WGS Group"))

lipo_comp = list(c("A", "B.1"), c("A", "B.2"), c("A", "C"), c("B.1", "B.2"), c("B.1", "C"), c("B.2", "C"))
p6.2b <- erp_stats %>% ggstripchart(x = "OspC_Type", y = "erp_num", color = "wgs_group") + 
  #facet_grid(~Region) +
  labs(x = "OspC Type", y = "Number of Erps") + 
  theme(legend.position = "none")

lipo_comp = list(c("A", "B.1"), c("A", "B.2"), c("A", "C"), c("B.1", "B.2"), c("B.1", "C"), c("B.2", "C"))
p6.3b <- mlp_stats %>% ggstripchart(x = "OspC_Type", y = "mlp_num", color = "wgs_group") + 
  #facet_grid(~Region) +
  labs(x = "OspC Type", y = "Number of Mlps") + 
  theme(legend.position = "none")

plot_grid(p6.1b, p6.2b, p6.3b, nrow = 3)
ggsave("../results/figures/paper/revision/Figure_6F.jpg", height = 6, width = 6)

# make correlation plots
plot_grid(p6.0c, p6.1c, p6.2c, p6.3c, nrow = 2)
ggsave("../results/figures/paper/revision/Figure_S7D.jpg", height = 6, width = 6)

lipo_stats %>% group_by(wgs_group) %>% summarize(dissem_prob = mean(dissem_bin, na.rm=T), orf_num = mean(lipo_num)) %>% 
  ggscatter(x = "orf_num", y = "dissem_prob")

p6.1d <- lipo_stats %>% group_by(OspC_Type) %>% summarize(dissem_prob = mean(dissem_bin, na.rm=T), orf_num = mean(lipo_num)) %>% 
  ggscatter(x = "orf_num", y = "dissem_prob", label = "OspC_Type", add = "reg.line", conf.int = TRUE) +
  stat_cor() + 
  labs(x = "Mean Number of Lipoproteins", y = "Mean Probability of Dissemination") + 
  coord_cartesian(ylim = c(0, 0.7))

p6.0d <- genome_stats %>% group_by(OspC_Type) %>% summarize(dissem_prob = mean(dissem_bin, na.rm=T), orf_num = mean(orf_num)) %>% 
  ggscatter(x = "orf_num", y = "dissem_prob", label = "OspC_Type", add = "reg.line", conf.int = TRUE) +
  stat_cor() + 
  labs(x = "Mean Number of ORFs", y = "Mean Probability of Dissemination") + 
  coord_cartesian(ylim = c(0, 0.7))

p6.2d <- erp_stats %>% group_by(OspC_Type) %>% summarize(dissem_prob = mean(dissem_bin, na.rm=T), orf_num = mean(erp_num)) %>% 
  ggscatter(x = "orf_num", y = "dissem_prob", label = "OspC_Type", add = "reg.line", conf.int = TRUE) +
  stat_cor() + 
  labs(x = "Mean Number of Erps", y = "Mean Probability of Dissemination") + 
  coord_cartesian(ylim = c(0, 0.7))

p6.3d <- mlp_stats %>% group_by(OspC_Type) %>% summarize(dissem_prob = mean(dissem_bin, na.rm=T), orf_num = mean(mlp_num)) %>% 
  ggscatter(x = "orf_num", y = "dissem_prob", label = "OspC_Type", add = "reg.line", conf.int = TRUE) +
  stat_cor() + 
  labs(x = "Mean Number of Mlps", y = "Mean Probability of Dissemination") + 
  coord_cartesian(ylim = c(0, 0.7))

plot_grid(p6.0d, p6.1d, p6.2d, p6.3d, nrow = 2)
ggsave("../results/figures/paper/revision/Figure_S7E.jpg", height = 6, width = 6)


# shade in OspC type A clade

# calculate correlation between different alleles

cor.test(as.numeric(gt.annotated[grep("g_1807", rownames(gt.annotated)),]), as.numeric(gt.annotated[grep("g_1021", rownames(gt.annotated)),])) 
cor.test(as.numeric(gt.annotated[grep("g_1820", rownames(gt.annotated)),]), as.numeric(gt.annotated[grep("g_1021", rownames(gt.annotated)),])) 

# Create plasmids PAM from PF-32 (Figure 4a) ####
plasmids <- data.frame(read_tsv("../plasmid_types_v3.tsv"))
rownames(plasmids) <- plasmids[,1]
plasmids <- plasmids[,-c(1)]

names(plasmids) <- sapply(names(plasmids), function(x) tail(strsplit(x, "_")[[1]], 1))
plasmids[is.na(plasmids)] <- 0
plasmids <- plasmids %>% select(-c(unknown))

# relabel tree to highlight clades
wgs_tree_HL <- wgs_tree + 
  geom_tippoint(aes(color=OspC_Type, subset=!is.na(OspC_Type)), size = 2) + 
  geom_highlight(node = 381, alpha = 0.5, fill = "darkgreen") + 
  geom_highlight(node = 379, fill = "steelblue", alpha = 0.2) 
#geom_cladelabel(node = 381, label = "")

OspC_type_tree_HL <- ggtree(x) %<+% metadata + 
  geom_tippoint(aes(color=OspC_Type, subset=!is.na(OspC_Type)), size = 2) + 
  geom_highlight(node = 381, alpha = 0.5, fill = "darkgreen") + 
  geom_highlight(node = 379, fill = "steelblue", alpha = 0.2) + 
  geom_treescale(x = 0.001, y = -25, fontsize = 6, linesize = 1, width = 0.002, offset = 2)
  #geom_cladelabel(node = 381, label = "")
OspC_type_tree_HL
p_plas <- heatMapPAM(location_tree,plasmids, col_colours="blue", colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)
ggsave("../results/figures/core_genome_plasmids_location.jpg", height = 12, width = 20)
# make sure you edit plasmids tsv file to replace | with _ so the columns read in correctly!
p_plas <- heatMapPAM(location_tree,plasmids, col_colours="blue", colnames_angle = -90,hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)
ggsave("../results/figures/core_genome_plasmids_location.jpg", height = 12, width = 20)

p_plas2 <- heatMapPAM(wgs_tree,plasmids, col_colours="blue", colnames_angle = -90,hjust =1, 
                      width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)
ggsave("../results/figures/core_genome_plasmids_wgs_group.jpg", height = 12, width = 20)

p_plas3 <- heatMapPAM(RST_type_tree,plasmids, col_colours="blue", colnames_angle = -90,
                      hjust =1, width=2, font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-50, border_colour=NULL)
ggsave("../results/figures/core_genome_plasmids_RST.jpg", height = 12, width = 20)

p_plas4 <- heatMapPAMpf32(OspC_type_tree_HL,plasmids, col_colours="blue", colnames_angle = -90,hjust =1, width=7, 
                          font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-40, border_colour=NULL)
ggsave("../results/figures/paper/revision/Figure_4A.jpg", height = 5, width = 8)


# Create plasmid PAM by alignment (Figure 4a_alt) #### 
plas2 <- data.frame(read_tsv("../plasmid_alignment_ref.tsv"))
rownames(plas2) <- plas2$Replicon
plas2 <- t(plas2[,-c(300:301)]) # replicon and length column removal is hard coded
colnames(plas2)[1] <- "Chrom"
p_plas5 <- heatMapPAMplas(OspC_type_tree_HL,plas2, col_colours="blue", colnames_angle = -90,hjust =1, width=7, 
                          font.size=4, cluster_cols=TRUE, null_colour="white",colnames_offset_y=-40, border_colour=NULL)
ggsave("../results/figures/paper/revision/Figure_S5.jpg", height = 5, width = 8)


# Calculate some associations ####
## Calculate associations between covariates
fisher.test(metadata_extended$WGS_Group, metadata_extended$Region) #, simulate.p.value = TRUE, B = 1e6)

#fisher.test(metadata_extended$Disseminated, metadata_extended$WGS_Group, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata_extended$RST_Type, metadata_extended$group, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata_extended$group, metadata_extended$OspC_Type, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata_extended$group, metadata_extended$Location, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata_extended$group, metadata_extended$Region, simulate.p.value = TRUE, B = 1e6)

# (no fig) plot of genome size ####
pGLength <- ggplot(metadata_extended, aes(x =groups, y = ORF_length)) + geom_jitter(width = 0.1) + 
                                                stat_summary(
                                                  fun.y = mean,
                                                  geom = "point",
                                                  size = 2,
                                                  color = "red"
                                                ) +
                                                stat_summary(fun.data = mean_se,
                                                             geom = "errorbar", color = "red",
                                                             width = 0.25) + 
                                                
                                                theme_bw() + 
                                                labs(x = "WGS Group", y = "Number of ORFs")
pGLength

# Figure 3 ####

pDissem <- ggplot(metadata_extended %>% filter(!is.na(Country)), aes(x =WGS_Group, y = dissem_bin, color = WGS_Group)) + geom_jitter(height = 0.02) + 
  stat_summary(
    fun.y = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  stat_summary(fun.y = function(x) BinomCI(sum(x), length(x))[1], 
               fun.ymin = function(x) BinomCI(sum(x), length(x))[2],
               fun.ymax = function(x) BinomCI(sum(x), length(x))[3],
               geom = "errorbar", color = "black",
               width = 0.25) + 
  theme_bw() + 
  labs(x = "WGS Group", y = "Probability of Dissemination") + 
  facet_grid(~Region) + 
  theme(legend.position = "none", axis.text=element_text(size = 18), axis.title=element_text(size = 18))
pDissem

dissem_tree <- ggtree(x) %<+% metadata_extended + 
  geom_tippoint(aes(color=Disseminated, subset=!is.na(Disseminated)), size = 2, alpha = 0.6) + 
  scale_color_discrete(labels = c("Disseminated", "Localized"), name = "") + 
  theme(legend.position = "left", legend.text=element_text(size=18))
#dissem_tree


my_comparisons = list( c("A","B"),c("A", "C"), c("B", "C"))
my_comparisons3 = list(c("A", "B"))
ORF_group <- metadata_extended %>% filter(!is.na(Country)) %>% 
  ggstripchart(x = "WGS_Group", y = "ORF_length", add="jitter", 
               color = "WGS_Group",
               add.params = list(size = 0.5, alpha =0.5)) + 
  #stat_compare_means(label = "p.signif", label.x=2) + 
  stat_compare_means(comparisons = my_comparisons3, label = "p.signif", tip.length = 0.01) + 
  labs(x = "WGS Group", y = "Number of ORFs") + 
  font("xlab", size = 20) + font("ylab", size = 20) + font("xy.text", size = 20)  + 
  facet_grid(~Region) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0.01) + 
  theme(legend.position = "none")
my_comparisons2 = list(c("US Midwest", "US Northeast"),c("US Midwest","EU Slovenia"), c("US Northeast", "EU Slovenia"))
my_comparisons4 = list(c("US Midwest", "US Northeast"))
ORF_region <- metadata_extended %>% 
  mutate(Country = recode(Region, Slovenia = "EU")) %>%filter(!is.na(Country)) %>% ggstripchart(x = "Country", y = "ORF_length", add="jitter", 
                                                                                                color = "WGS_Group",
                                                                                                add.params = list(size = 0.5, alpha =0.5)) + 
  #stat_compare_means(label = "p.signif", label.x=2) + 
  labs(x = "", y = "Number of ORFs") + 
  font("xlab", size = 20) + font("ylab", size = 20) + font("xy.text", size = 20)  + 
  facet_grid(~WGS_Group) + 
  stat_compare_means(comparisons = my_comparisons2, label = "p.signif") + 
  stat_compare_means(comparisons = my_comparisons4, label = "p.signif") + 
  theme(axis.text.x = element_text(angle = -90, size = 12), legend.position = "none")

#ORF_region
plot_grid( ORF_region,ORF_group, pDissem, labels = c("A", "B", "C"), nrow = 1)#, rel_widths = c(2,1.2,2))
ggsave("../results/figures/paper/revision/Figure_3.jpg", height = 5, width = 12)

# Figure S3A ####
pDissem2 <- ggplot(metadata_extended %>% filter(!is.na(Country)), aes(x =wgs_group, y = dissem_bin, color = wgs_group)) + geom_jitter(height = 0.02) + 
  stat_summary(
    fun.y = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  stat_summary(fun.y = function(x) BinomCI(sum(x), length(x))[1], 
               fun.ymin = function(x) BinomCI(sum(x), length(x))[2],
               fun.ymax = function(x) BinomCI(sum(x), length(x))[3],
               geom = "errorbar", color = "black",
               width = 0.25) + 
  theme_bw() + 
  labs(x = "WGS Group", y = "Probability of Dissemination") + 
  facet_grid(~Country) + 
  theme(legend.position = "none", axis.text=element_text(size = 18), axis.title=element_text(size = 18))



pDissem2.1 <- ggplot(metadata_extended %>% filter(!is.na(Country)), aes(x =WGS_Group, y = dissem_bin, color = wgs_group)) + geom_jitter(height = 0.02) + 
  stat_summary(
    fun.y = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  stat_summary(fun.y = function(x) BinomCI(sum(x), length(x))[1], 
               fun.ymin = function(x) BinomCI(sum(x), length(x))[2],
               fun.ymax = function(x) BinomCI(sum(x), length(x))[3],
               geom = "errorbar", color = "black",
               width = 0.25) + 
  theme_bw() + 
  labs(x = "WGS Group", y = "Probability of Dissemination") + 
  facet_grid(~Region) + 
  theme(legend.position = "none", axis.text=element_text(size = 18), axis.title=element_text(size = 18))

pDissem2.2 <- ggplot(metadata_extended %>% filter(!is.na(Country)), aes(x =wgs_group, y = dissem_bin, color = wgs_group)) + geom_jitter(height = 0.02) + 
  stat_summary(
    fun.y = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  stat_summary(fun.y = function(x) BinomCI(sum(x), length(x))[1], 
               fun.ymin = function(x) BinomCI(sum(x), length(x))[2],
               fun.ymax = function(x) BinomCI(sum(x), length(x))[3],
               geom = "errorbar", color = "black",
               width = 0.25) + 
  theme_bw() + 
  labs(x = "WGS Group", y = "Probability of Dissemination") + 
  facet_grid(~Region) + 
  theme(legend.position = "none", axis.text=element_text(size = 18), axis.title=element_text(size = 18))

my_comparisons = list(c("a","c"))#, c("B","C"), c("A", "C"))
my_comparisons3 = list(c("a","c"))
ORF_group <- metadata_extended %>% filter(!is.na(Country)) %>% ggstripchart(x = "wgs_group", y = "ORF_length", add="jitter", 
                                                                            color = "wgs_group",
                                                                            add.params = list(size = 0.5, alpha =0.5)) + 
  #stat_compare_means(label = "p.signif", label.x=2) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(comparisons = my_comparisons3, label = "p.signif") + 
  labs(x = "WGS Group", y = "Number of ORFs") + 
  font("xlab", size = 20) + font("ylab", size = 20) + font("xy.text", size = 20)  + 
  facet_grid(~Country)
my_comparisons2 = list(c("US","Slovenia"))
ORF_region <- metadata_extended %>% filter(!is.na(Country)) %>% ggstripchart(x = "Country", y = "ORF_length", add="jitter", 
                                                                             color = "wgs_group",
                                                                             add.params = list(size = 0.5, alpha =0.5)) + 
  #stat_compare_means(label = "p.signif", label.x=2) + 
  stat_compare_means(comparisons = my_comparisons2, label = "p.signif") + 
  labs(x = "", y = "Number of ORFs") + 
  font("xlab", size = 20) + font("ylab", size = 20) + font("xy.text", size = 20)  + 
  facet_grid(~wgs_group)
#ORF_region
wgs2_tree_HL <- ggtree(x) %<+% metadata_extended + 
  geom_tippoint(aes(color=wgs_group, subset=!is.na(wgs_group)), size = 4) + 
  geom_highlight(node = 381, alpha = 0.5, fill = "darkgreen") + 
  geom_highlight(node = 379, fill = "steelblue", alpha = 0.2) 
  #geom_cladelabel(node = 381, label = "")


p_hm2b <- heatMapPAM(wgs2_tree_HL,t(genotype.trunc2), colnames=FALSE, col_colours="blue", 
                     colnames_angle = 90,hjust =0, width=10, font.size=5, cluster_cols=TRUE, 
                     null_colour="white",colnames_offset_y=-20, border_colour=NULL,
                     colnames_position = "bottom")


ggsave("../results/figures/paper/revision/Figure_S3A.jpg", height = 5, width = 10)

# calculate dissemination rates:
metadata_extended %>% group_by(wgs_group, Country) %>% summarize(dissem = mean(dissem_bin, na.rm=T))
metadata_extended %>% group_by(WGS_Group, Country) %>% summarize(dissem = mean(dissem_bin, na.rm=T))

mt_s <- metadata_extended %>% filter(Country == "US")

# calculate dissemination statistics
fisher.test(mt_s$Disseminated, mt_s$WGS_Group, simulate.p.value = TRUE, B = 1e6)
fisher.test(mt_s$Disseminated, mt_s$wgs_group, simulate.p.value = TRUE, B = 1e6)

fisher.test(metadata_extended$Disseminated, metadata_extended$WGS_Group, simulate.p.value = TRUE, B = 1e6)
fisher.test(metadata_extended$Disseminated, metadata_extended$wgs_group, simulate.p.value = TRUE, B = 1e6)

# Slovenian vs US
fisher.test(metadata_extended$Disseminated, metadata_extended$Country, simulate.p.value = TRUE, B = 1e6)

metadata_extended %>% group_by(Country) %>% summarize(dissem_rate = mean(dissem_bin, na.rm=T))

metadata_extended %>% group_by(Country, WGS_Group) %>% summarize(dissem_rate = mean(dissem_bin, na.rm=T))
metadata_extended %>% group_by(Country, wgs_group) %>% summarize(dissem_rate = mean(dissem_bin, na.rm=T))

# RST type vs dissem
fisher.test(metadata_extended$Disseminated, metadata_extended$RST_Type, simulate.p.value = TRUE, B = 1e6)
fisher.test(metadata_extended$Disseminated, metadata_extended$OspC_Type, simulate.p.value = TRUE, B = 1e6)
fisher.test(metadata_extended$Disseminated, metadata_extended$OspCA, simulate.p.value = TRUE, B = 1e6)

chisq.test(metadata_extended$Disseminated, metadata_extended$RST_Type, simulate.p.value = TRUE, B = 1e6)
chisq.test(metadata_extended$Disseminated, metadata_extended$OspC_Type, simulate.p.value = TRUE, B = 1e6)
chisq.test(metadata_extended$Disseminated, metadata_extended$OspCA, simulate.p.value = TRUE, B = 1e6)


#### Construct contigency tables

# Slovenian vs US
  
tableDissem <- function(col1, col2, rowname_prefix){
  outdf <- as.data.frame.matrix(table(col1, col2))
  outdf$fisher_test <- c(round(fisher.test(col1,col2, simulate.p.value = TRUE, B = 1e6)$p.value, 6), rep(NA, nrow(outdf)-1))
  outdf$chisq_test <- c(round(chisq.test(col1,col2)$p.value, 6), rep(NA, nrow(outdf)-1))
  rownames(outdf) <- paste(rep(rowname_prefix, nrow(outdf)), rownames(outdf))
  outdf
}

dissem_stats <- tableDissem(mt_s$WGS_Group, mt_s$Disseminated, "(US isolates only) WGS Group")
dissem_stats <- rbind(dissem_stats, tableDissem(mt_s$wgs_group, mt_s$Disseminated, "(US isolates only) WGS Group (split B)"))
dissem_stats <- rbind(dissem_stats, tableDissem(metadata_extended$WGS_Group,metadata_extended$Disseminated, "WGS Group"))
dissem_stats <- rbind(dissem_stats, tableDissem(metadata_extended$wgs_group,metadata_extended$Disseminated, "WGS Group (split B)"))
dissem_stats <- rbind(dissem_stats, tableDissem(metadata_extended$Country, metadata_extended$Disseminated, "Country"))
dissem_stats <- rbind(dissem_stats, tableDissem(metadata_extended$RST_Type, metadata_extended$Disseminated, "RST Type"))
dissem_stats <- rbind(dissem_stats, tableDissem(metadata_extended$OspC_Type, metadata_extended$Disseminated, "Osp C Type"))
dissem_stats <- rbind(dissem_stats, tableDissem(metadata_extended$OspCA, metadata_extended$Disseminated, "Osp C Type A binary"))
colnames(dissem_stats) <- c("Number of Disseminated Isolates", "Number of Localized Isolates", "Fisher Exact test P value", "Chi Square test P value")
write.csv(dissem_stats, "../results/tables/revision/ST_4_new.csv", na = "")

# Figure S3B-D ####
plot_grid(mds_plot2, pDissem2, labels = c("B", "C"),nrow=1)
ggsave("../results/figures/paper/revision/Figure_3B-C.jpg",height = 5, width = 12)

fisher.test(metadata_extended$Severity_Score, metadata_extended$WGS_Group, simulate.p.value = TRUE, B = 1e6)

# Test dissem US vs Slovenia ####
pDissem3 <- ggplot(metadata_extended %>% filter(!is.na(Country)), aes(x =Country, y = dissem_bin)) + geom_jitter(height = 0.02) + 
  stat_summary(
    fun.y = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  stat_summary(fun.y = function(x) BinomCI(sum(x), length(x))[1], 
               fun.ymin = function(x) BinomCI(sum(x), length(x))[2],
               fun.ymax = function(x) BinomCI(sum(x), length(x))[3],
               geom = "errorbar", color = "black",
               width = 0.25) + 
  theme_bw() + 
  labs(x = "WGS Group", y = "Probability of Dissemination") + 
  theme(legend.position = "none", axis.text=element_text(size = 18), axis.title=element_text(size = 18))
pDissem3
fisher.test(metadata_extended$Disseminated, metadata_extended$Country, simulate.p.value = TRUE, B = 1e6)

# Figure S4 ####
pDissem4 <- ggplot(metadata_extended %>% filter(!is.na(OspC_Type)), aes(x =OspC_Type, y = dissem_bin, color=wgs_group)) + geom_jitter(height = 0.02) + 
  stat_summary(
    fun.y = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  stat_summary(fun.y = function(x) BinomCI(sum(x), length(x))[1], 
               fun.ymin = function(x) BinomCI(sum(x), length(x))[2],
               fun.ymax = function(x) BinomCI(sum(x), length(x))[3],
               geom = "errorbar", color = "black",
               width = 0.25) + 
  theme_bw() + 
  labs(x = "OspC Type", y = "Probability of Dissemination") + 
  theme(legend.position = "top", axis.text=element_text(size = 18), axis.title=element_text(size = 18), 
        legend.text = element_text(size = 24), legend.title = element_text(size = 24)) + 
  facet_grid(rows=vars(Country)) + 
  guides(color = guide_legend(override.aes = list(size=5)))
fisher.test(metadata_extended$Disseminated, metadata_extended$OspC_Type, simulate.p.value = TRUE, B = 1e6)

pDissem5 <- ggplot(metadata_extended %>% filter(!is.na(OspC_Type)), aes(x =RST_Type, y = dissem_bin, color=wgs_group)) + geom_jitter(height = 0.02) + 
  stat_summary(
    fun.y = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  stat_summary(fun.y = function(x) BinomCI(sum(x), length(x))[1], 
               fun.ymin = function(x) BinomCI(sum(x), length(x))[2],
               fun.ymax = function(x) BinomCI(sum(x), length(x))[3],
               geom = "errorbar", color = "black",
               width = 0.25) + 
  theme_bw() + 
  labs(x = "RST", y = "Probability of Dissemination") + 
  theme(legend.position = "none", axis.text=element_text(size = 18), axis.title=element_text(size = 18)) +
  facet_grid(~Country)

#mycomps <- list(c("1", "2"), c("1", "3"), c("2", "3"))
pSevere <- ggplot(metadata_extended, aes(x =wgs_group, y = Severity_Score, color = wgs_group)) + geom_jitter(height = 0.02) + 
  stat_summary(
    fun.y = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar", color = "black",
               width = 0.25) + 
  #stat_compare_means(comparisons = mycomps, tip.length = 0.01) + 
  theme_bw() + 
  #facet_grid(~Region) +
  labs(x = "WGS Group", y = "Severity") + 
  theme(legend.position = "none", axis.text=element_text(size = 18), axis.title=element_text(size = 18))
pSev <- metadata_extended %>% 
  filter(!is.na(Severity)) %>% 
  ggplot(aes(x = OspC_Type, fill = wgs_group)) + geom_bar() + facet_grid(rows = vars(Severity)) +
  theme_bw() +
  theme(legend.position = "none", axis.text=element_text(size = 18), axis.title=element_text(size = 18))


pg1 <- plot_grid(pDissem4, labels = c("A"))
pg2 <- plot_grid(pDissem5, pSev, labels = c("B", "C"), rel_widths = c(1,1))
plot_grid(pg1, pg2, nrow = 2)
fisher.test(metadata_extended$Disseminated, metadata_extended$RST_Type, simulate.p.value = TRUE, B = 1e6)
ggsave("../results/figures/paper/revision/Figure_S4.jpg",height = 8, width = 8)


# Fig S4B accessory genome PAM heatmap #### 

accessory_tree <- ggtree(a) %<+% metadata_extended + 
  geom_tippoint(aes(color=wgs_group, subset=!is.na(wgs_group)), size = 2)
accessory_tree

p_hm_ac <- heatMapPAM(accessory_tree,t(genotype.trunc2), colnames=FALSE, col_colours="blue", 
                      colnames_angle = 90,hjust =0, width=10, font.size=5, cluster_cols=TRUE, 
                      null_colour="white",colnames_offset_y=-20, border_colour=NULL,
                      colnames_position = "bottom", presence_label = "ORF Present", absence_label = "ORF Absent") 
ggsave("../results/figures/paper/revision/Figure_S4A.jpg", height = 5, width = 6)

p_hm_ac_p <- heatMapPAM(accessory_tree,plasmids,, col_colours="blue", colnames_angle = -90,
                        hjust =1, width=10, font.size=2, cluster_cols=TRUE, 
                        null_colour="white",colnames_offset_y=-65, border_colour=NULL,
                        presence_label = "PF-32 Present", absence_label = "PF-32 Absent")
ggsave("../results/figures/paper/revision/Figure_S6B.jpg", height = 5, width = 6)

# some modelling scrap #####

## my_comparisons4 = list(c("D","L"))
metadata_extended %>% filter(!is.na(Disseminated)) %>% ggstripchart(x = "Disseminated", y = "ORF_length", add="jitter", 
                                                                    add.params = list(size = 0.5, alpha =0.5)) + 
  #stat_compare_means(label = "p.signif", label.x=2) + 
  stat_compare_means(comparisons = my_comparisons4, label = "p.signif") + 
  labs(x = "", y = "Number of ORFs") + 
  font("xlab", size = 20) + font("ylab", size = 20) + font("xy.text", size = 20)  + 
  facet_grid(~Country~wgs_group)

## 
# analysis by fine-grain geography
NE_samples <- metadata_extended %>% filter(Region == "US Northeast")

fisher.test(NE_samples$WGS_Group, NE_samples$Location, simulate.p.value = TRUE, B = 1e6)

# Figure 1 ####

# calculate associations between metadata variables
#ggpairs(metadata_extended %>% filter(c(OspC_Type, RST_Type, groups)), cardinality_threshold = 25)
#ggsave("../results/figures/OspC_Type_vs_RST_pairs.jpg")
pA <- metadata %>% 
  filter(!is.na(RST_Type)) %>% 
  ggplot(aes(x = OspC_Type, fill=RST_Type)) + geom_bar() + facet_grid(rows = vars(Region)) +
  theme_bw() +
  theme(legend.position="top")+
  scale_fill_manual(values = c("blue", "orange", "green"))

pB <- metadata_extended %>% 
  filter(!is.na(RST_Type)) %>% 
  ggplot(aes(x = OspC_Type, fill=WGS_Group)) + geom_bar() + facet_grid(rows = vars(Region)) + 
  theme_bw() + 
  theme(legend.position = "top") 
#plot_grid(pA, pB, labels = c("A", "B"))

# Plot and color by groups
mds_plot <- ggscatter(mds, x = "V1", y = "V2", 
                      #label = rownames(d),
                      color = "WGS_Group",
                      size = 1, 
                      ellipse = TRUE,
                      #ellipse.type = "convex",
                      repel = TRUE) + labs(x = "MDS Coordinate 1", y = "MDS Coordinate 2")

mds_RST <- ggscatter(metadata_extended %>% filter(!is.na(RST_Type)), x = "V1", y = "V2", 
                     #label = rownames(d),
                     color = "RST_Type",
                     size = 1,
                     palette = c("blue", "orange", "green"),
                     ellipse=TRUE,
                     repel = TRUE) + labs(x = "MDS Coordinate 1", y = "MDS Coordinate 2")

plot_grid(pA, pB, mds_RST, mds_plot, labels = c("A", "B", "C", "D"), nrow = 2)

ggsave("../results/figures/paper/revision/Figure_1.jpg", height = 10, width = 10)
#ggsave("../results/figures/OspC_Type_by_Region.jpg")

# More associations scrap ####

#metadata %>% 
#    filter(!is.na(RST_Type)) %>% 
#    ggplot(aes(x = Disseminated)) + geom_bar() + facet_grid(rows = vars(RST_Type))
#ggsave("../results/figures/Disseminated_vs_RST_Type.jpg")
#metadata %>% 
#filter(!is.na(Disseminated)) %>% 
#ggplot(aes(x = OspC_Type)) + geom_bar() + facet_grid(rows = vars(Disseminated))
#ggsave("../results/figures/OspC_Type_vs_Dissemination.jpg")

#chisq.test(metadata$Location, metadata$OspC_Type)
#chisq.test(metadata$Disseminated, metadata$RST_Type)
#chisq.test(metadata$Disseminated, metadata$OspC_Type)
#chisq.test(metadata$Severity, metadata$RST_Type)
#chisq.test(metadata$Severity, metadata$OspC_Type)
#chisq.test(metadata$Disseminated, metadata$OspCA)
#chisq.test(metadata$RST_Type, metadata$OspC_Type)
#chisq.test(metadata$Location, metadata$OspC_Type)


## Calculate associations between categorical covariates

#fisher.test(metadata_extended$Disseminated, metadata_extended$RST_Type, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata_extended$WGS_Group, metadata_extended$Region)#, simulate.p.value = TRUE, B = 1e6)
fisher.test(metadata_extended$RST_Type, metadata_extended$Region)#, simulate.p.value = TRUE, B = 1e6)
fisher.test(metadata_extended$OspC_Type, metadata_extended$Region, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata$Disseminated, metadata$OspC_Type, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata$Severity, metadata$RST_Type, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata$Severity, metadata$OspC_Type, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata$Disseminated, metadata$OspCA, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata$RST_Type, metadata$OspC_Type, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata$Location, metadata$OspC_Type, simulate.p.value = TRUE, B = 1e6)
#fisher.test(metadata$Location, metadata$RST_Type, simulate.p.value = TRUE, B = 1e6)

US_isolates <- metadata_extended %>% filter(Country == "US")
EU_isolates <- metadata_extended %>% filter(Country == "Slovenia")

fisher.test(US_isolates$Disseminated, US_isolates$WGS_Group, simulate.p.value = TRUE, B = 1e6)
fisher.test(EU_isolates$Disseminated, EU_isolates$WGS_Group, simulate.p.value = TRUE, B = 1e6)

## 

# Figure 2 ####

# generate face-to-face tree of clonalframeML tree
p1 <- ggtree(x2) %<+% metadata_extended
p2 <- ggtree(x1) %<+% metadata_extended
d1 <- p1$data
d2 <- p2$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
pp <- p1 + geom_tree(data=d2)
d1 <- p1$data

dd <- bind_rows(d1, d2) %>% filter(isTip == TRUE)

f2fa <- pp + geom_line(aes(x, y, group=label, color=wgs_group), data=dd, alpha = 0.6, size = 2) #+ 
#geom_tippoint(aes(color=wgs_group)) 
#geom_tippoint(data=dd, aes(color=wgs_group))
#f2fa
plot_grid(p1, p2)


### Generate figure of face-to-face core and accessory tress (currently S4B)

p1 <- ggtree(x) %<+% metadata_extended
p2 <- ggtree(a) %<+% metadata_extended

d1 <- p1$data
d2 <- p2$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
p1$data$x <- p1$data$x*30
pp <- p1 + geom_tree(data=d2)
d1 <- p1$data

dd <- bind_rows(d1, d2) %>% filter(isTip == TRUE)

f2f <- pp + geom_line(aes(x, y, group=label, color=wgs_group), data=dd, alpha = 0.6, size = 2) #+ 
#geom_tippoint(aes(color=wgs_group)) 
#geom_tippoint(data=dd, aes(color=wgs_group))
f2f
#ggsave("../results/figures/paper/revision/Figure_S4B_core_and_accessory_genome_trees_face_to_face.jpg", height = 4, width = 4)

# Figure S2C (ML vs MCC F2F) ####

p1 <- ggtree(b) %<+% (metadata_extended %>% rename(`WGS Group` = wgs_group )) + 
  geom_treescale(y = -10, offset = 2)
p2 <- ggtree(iq_tree) %<+% (metadata_extended %>% rename(`WGS Group` = wgs_group ))
d1 <- p1$data
d2 <- p2$data
## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
expansion_factor <- 1
d2$x <- max(d2$x*expansion_factor) - d2$x*expansion_factor + max(d1$x*expansion_factor)
p1$data$x <- p1$data$x*expansion_factor
pp <- p1 + geom_tree(data=d2)
d1 <- p1$data

dd <- bind_rows(d1, d2) %>% filter(!is.na(OspC_Type)) %>% filter(isTip == TRUE)
f2f2 <- pp + geom_line(aes(x, y, group=label, color=`WGS Group`), data=dd, alpha = 0.6, size = 2) + 
  new_scale_color() + 
  geom_nodepoint(aes(color = posterior, subset = posterior > 0.9)) + 
  new_scale_color()+ 
  scale_colour_gradient2() +
  geom_nodepoint(data = d2, aes(color = UFboot, subset = UFboot > 90))
f2f2 + ggtree::vexpand(0.2, 1)
f2f2
ggsave("../results/figures/paper/revision/Figure_S2C.jpg", height = 8, width = 8)

# Figure S2D (MCC vs OspC F2F) ####

p1 <- ggtree(b) %<+% (metadata_extended %>% rename(`WGS Group` = wgs_group )) + 
  geom_treescale(y = -10, offset = -30, width = 0.05, label = as.character(0.05 / 50))
p2 <- ggtree(O_tree) %<+% (metadata_extended %>% rename(`WGS Group` = wgs_group ))
d1 <- p1$data
d2 <- p2$data
## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
expansion_factor <- 50
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
p1$data$x <- p1$data$x*expansion_factor
pp <- p1 + geom_tree(data=d2)
d1 <- p1$data

dd <- bind_rows(d1, d2) %>% filter(!is.na(OspC_Type)) %>% filter(isTip == TRUE)
f2f3 <- pp + geom_line(aes(x, y, group=label, color=`OspC_Type`), data=dd, alpha = 0.6, size = 2) + 
  new_scale_color() + 
  geom_nodepoint(aes(color = posterior, subset = posterior > 0.9)) + 
  geom_nodepoint(data = d2, aes(color = posterior, subset = posterior > 0.9))
f2f3 + ggtree::vexpand(0.2, 1)
f2f3
ggsave("../results/figures/paper/revision/Figure_S2D.jpg", height = 8, width = 8)


##
## Figure 2

# make wgs tree with additional annotation
# make annotated accessory trees


acc_tree <- ggtree(a) %<+% metadata_extended + geom_tippoint(aes(color=WGS_Group, subset=!is.na(WGS_Group)), size = 4)
acc2_tree <- ggtree(a) %<+% metadata_extended + geom_tippoint(aes(color=wgs_group, subset=!is.na(wgs_group)), size = 4)
acc_tree_ospC <- ggtree(a) %<+% metadata_extended + geom_tippoint(aes(color=OspC_Type, subset=!is.na(OspC_Type)), size = 4)

# Make OspC tree annotated with WGS groups
ospC_tree_WGS <- ggtree(y) %<+% metadata_extended + 
  geom_tippoint(aes(color=WGS_Group, subset=!is.na(WGS_Group)), size = 4)


# assemble plot
#plot_grid(wgs_tree, wgs2_tree, location_tree, f2f, wgs_tree_ospc, ospC_tree, ospC_tree_WGS, f2f2, nrow = 2, labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
#ggsave("../results/figures/paper/revision/Figure_2.jpg", height = 8, width = 10)


# try to reduce the number of trees, as requested by the reviewers.
P <- ggtree(x, layout = "radial")

metadata_heatmap <- metadata_extended[c("label","wgs_group", "OspC_Type", "Region")]
rownames(metadata_heatmap) <- metadata_heatmap$label
metadata_heatmap <- subset(metadata_heatmap, select= -c(label))
heat1 <- subset(metadata_heatmap, select = c(wgs_group, Region))
colnames(heat1) <- c("WGS", "Region")
heat2 <- subset(metadata_heatmap, select = c(OspC_Type))
colnames(heat2) <- c("OspC")
h1 <- gheatmap(p, heat1, offset = 0.0035,
               legend_title = "WGS Group\n\nRegion",
               colnames_angle = 0, 
               colnames_offset_y = -4)
h2 <- h1 + new_scale_fill() 
h3 <- gheatmap(h2, heat2, width = 0.5, 
         legend_title = "OspC Type",
         colnames_angle = 0, 
         colnames_offset_y = -4) + 
  ggtree::vexpand(.1, -1)

# assemble new version of plot
# Supplementary figure 2E_F #### 
plot_grid(f2f, f2f2, nrow = 2, labels = c("E", "F"))
ggsave("../results/figures/paper/revision/Figure_S2_EF.jpg", height = 10, width = 7)

# Figure 2 Revised ####

md <- metadata_extended[,c("label", "Region", "RST_Type")] %>% 
  pivot_longer(!label)
md$pos <- paste(md$name, md$value)
md$pos <- gsub("Region ", "", md$pos)
md$pos <- gsub("_Type", "", md$pos)
md$pos <- factor(md$pos, levels = c("RST 1", "RST 2", "RST 3", "EU Slovenia", "US Midwest", "US Northeast"))
md <- md %>% mutate(pos2 = as.numeric(pos))

md2 <- metadata_extended[,c("label", "Disseminated")] %>% 
  filter(!is.na(Disseminated)) %>% 
  mutate(Disseminated = recode(Disseminated, D = "Disseminated", L = "Localized"))
names(md2)[2] <- "Clinical Status"

md3 <- metadata_extended[,c("label", "ORF_length", "wgs_group")]
names(md3)[3] <- "WGS Group"

q1 <- ggtree(b, layout="fan", open.angle=5) %<+% 
  metadata_extended[,c("label", "OspC_Type", "wgs_group")] + 
  geom_tiplab(aes(label = paste("OspC Type ",OspC_Type, sep="")), size=1.2, offset = 0.00015) +
  geom_nodelab(aes(label = round(posterior,2), subset = posterior > 0.9), 
               ,hjust = 1.3,alpha = 0.4, color = "blue") + 
  #geom_tippoint(aes(color = wgs_group), alpha = 0.4)+
  geom_treescale(x = 0.0015, y = 280, offset = 4, width = 0.001)+
  scale_fill_manual(values=hue_pal()(23), guide = "none") +
  geom_fruit(geom = geom_bar,
             mapping = aes(y=label,x=1,fill=OspC_Type),
             stat='identity',
             alpha = 0.5,
             pwidth = 0.18,
             offset = 0.03) +
  new_scale_fill()+ 
  geom_fruit(data = md, geom = geom_tile,
             mapping = aes(y=label,x=name,fill=pos),
             stat='identity',
             alpha = 0.5,
             offset = 0.02, pwidth = 0.06) + 
  scale_fill_viridis_d(option="D", name="Variable") +
  new_scale_fill() +
  new_scale_color() + 
  geom_fruit(data = md2, geom = geom_star,
             mapping = aes(y = label, fill = `Clinical Status`, starshape = `Clinical Status`,color = `Clinical Status`),
             stat='identity',
             size = 0.8,
             offset = 0.08) + 
  new_scale_color()+
  geom_fruit(geom = geom_point,
             data = md3,
             mapping = aes(y=label,x=1, color = `WGS Group`),
             stat='identity',
             pwidth = 0.05, 
             offset = 0,
             size = 0.8) 
  
  

q1
ggsave("../results/figures/paper/revision/Figure_2.jpg", height = 8, width = 8)

# Figure S2A (ML tree) ####
q3 <- ggtree(iq_tree, layout="fan", open.angle=5) %<+% 
  metadata_extended[,c("label", "OspC_Type", "wgs_group")] + 
  geom_tiplab(aes(label = paste("OspC Type ",OspC_Type, sep=""), color = OspC_Type), size=1.2, offset = 0.00015) +
  geom_nodelab(aes(label = round(UFboot,2), subset = UFboot > 90), 
               ,hjust = 1.3,alpha = 0.4, color = "blue") + 
  #geom_tippoint(aes(color = wgs_group), alpha = 0.4)+
  geom_treescale(x = 0.0015, y = 280, offset = 4, width = 0.001)+
  scale_fill_manual(values=hue_pal()(23), guide = "none") +
  geom_fruit(geom = geom_bar,
             mapping = aes(y=label,x=1,fill=OspC_Type),
             stat='identity',
             alpha = 0.5,
             pwidth = 0.18,
             offset = 0.03) +
  new_scale_fill()+ 
  geom_fruit(data = md, geom = geom_tile,
             mapping = aes(y=label,x=name,fill=pos),
             stat='identity',
             alpha = 0.5,
             offset = 0.02, pwidth = 0.06) + 
  scale_fill_viridis_d(option="D", name="Variable") +
  new_scale_fill() +
  new_scale_color() + 
  geom_fruit(data = md2, geom = geom_star,
             mapping = aes(y = label, fill = `Clinical Status`, starshape = `Clinical Status`,color = `Clinical Status`),
             stat='identity',
             size = 0.8,
             offset = 0.08) + 
  new_scale_color()+
  geom_fruit(geom = geom_point,
             data = md3,
             mapping = aes(y=label,x=1, color = `WGS Group`),
             stat='identity',
             pwidth = 0.05, 
             offset = 0,
             size = 0.8) 



q3
ggsave("../results/figures/paper/revision/Figure_S2A.jpg", height = 8, width = 8)

# Figure S2E (OspC circular tree) ####
metadata_extended_ospC <- metadata_extended[,c("label", "OspC_Type", "wgs_group")]
metadata_extended_ospC <- rbind(metadata_extended_ospC, 
                                c("URI20_1", "G", "B.2"),
                                c("URI20_2", "G", "B.2"))

q2 <- ggtree(O_tree, layout="fan", open.angle=5) %<+% 
  metadata_extended_ospC + 
  geom_tiplab(aes(label = paste("OspC Type ",OspC_Type, sep="")), size=1.2, offset = 0.0075) +
  geom_nodelab(aes(label = round(posterior,2), subset = posterior > 0.9), 
               ,hjust = 1.3,alpha = 0.4, color = "blue") + 
  #geom_tippoint(aes(color = wgs_group), alpha = 0.4)+
  scale_fill_manual(values=hue_pal()(23), guide = "none") +
  geom_fruit(geom = geom_bar,
             mapping = aes(y=label,x=1,fill=OspC_Type),
             stat='identity',
             alpha = 0.5,
             pwidth = 0.18,
             offset = 0.03) +
  new_scale_fill()+ 
  geom_fruit(data = md, geom = geom_tile,
             mapping = aes(y=label,x=name,fill=pos),
             stat='identity',
             alpha = 0.5,
             offset = 0.02, pwidth = 0.06) + 
  scale_fill_viridis_d(option="D", name="Variable") +
  new_scale_fill() +
  new_scale_color() + 
  geom_fruit(data = md2, geom = geom_star,
             mapping = aes(y = label, fill = `Clinical Status`, starshape = `Clinical Status`,color = `Clinical Status`),
             stat='identity',
             size = 0.8,
             offset = 0.08) + 
  new_scale_color()+
  geom_fruit(geom = geom_point,
             data = md3,
             mapping = aes(y=label,x=1, color = `WGS Group`),
             stat='identity',
             pwidth = 0.05, 
             offset = 0,
             size = 0.8) + 
  geom_treescale(x = 0.05, y = 280, offset = 4, width = 0.02)
q2
ggsave("../results/figures/paper/revision/Figure_S2E.jpg", height = 8, width = 8)


# temp figure ####

acc_tree_RST <- ggtree(a) %<+% metadata_extended + 
  geom_tippoint(aes(color=OspC_Type))
acc_tree_RST


names(metadata_extended)

metadata_extended %>% filter(OspC_Type == "S") %>% select(label, RST_Type, OspC_Type, wgs_group, MLST)

wgs_tree + geom_tiplab(aes(label=label)) + geom_nodelab(aes(label = node)) + 
  geom_cladelabel(node = 381, color = "green", label = "   OspC Type A") +
  geom_cladelabel(node = 379, color = "red", label = "RST1") + 
  geom_highlight(node = 381, fill = "orange") + 
  geom_highlight(node = 379)
  
ggsave("../results/figures/paper/revision/Figure_annotated.jpg", height = 50, width = 8, limitsize=FALSE)


# Fig 4C-D PF32 presence/absence associations with dissemination ####
# this analysis uses the Pfam32 gene analysis

# create merged dataset
plasmid_merge <- plasmids
plasmid_merge <- plasmid_merge[, colSums(plasmid_merge) < nrow(plasmid_merge)]
nplasmids <- ncol(plasmid_merge)
plasmid_merge$label <- rownames(plasmids)
plasmid_merge <- inner_join(plasmid_merge, metadata_extended)

# need to deal with the plasmids being hard-coded...

# conduct association tests

plasnames <- rep(NA, nplasmids)
pvals <- rep(NA, nplasmids)
estim <- rep(NA, nplasmids)
lowerint <- rep(NA, nplasmids)
upperint <- rep(NA, nplasmids)
for(i in 1:length(plasnames)){
  test_tmp <- fisher.test(plasmid_merge[,i], plasmid_merge$dissem_bin)
  plasnames[i] <- colnames(plasmid_merge)[i]
  pvals[i] <- test_tmp$p.value
  estim[i] <- test_tmp$estimate
  lowerint[i] <- test_tmp$conf.int[1]
  upperint[i] <- test_tmp$conf.int[2]
  #print(paste(i, test_tmp$p.value, " ", colnames(plasmid_merge)[i]))   
}
plasmid_assoc <- data.frame(plasmid = plasnames, pval = pvals, OR = estim, conflower = lowerint, confupper = upperint, logp = -1*log10(pvals),
                            adjustedP = p.adjust(pvals, method="fdr"))

# repeat for country association tests
plasmid_merge$country_bin <- ifelse(plasmid_merge$Country == "Slovenia", 1, 0)

plasnamesC <- rep(NA, nplasmids)
pvalsC <- rep(NA, nplasmids)
estimC <- rep(NA, nplasmids)
lowerintC <- rep(NA, nplasmids)
upperintC <- rep(NA, nplasmids)
for(i in 1:length(plasnamesC)){
  test_tmpC <- fisher.test(plasmid_merge[,i], plasmid_merge$country_bin)
  plasnamesC[i] <- colnames(plasmid_merge)[i]
  pvalsC[i] <- test_tmpC$p.value
  estimC[i] <- test_tmpC$estimate
  lowerintC[i] <- test_tmpC$conf.int[1]
  upperintC[i] <- test_tmpC$conf.int[2]
  #print(paste(i, test_tmp$p.value, " ", colnames(plasmid_merge)[i]))   
}
plasmid_assocC <- data.frame(plasmid = plasnamesC, pval = pvalsC, OR = estimC, conflower = lowerintC, confupper = upperintC, 
                             logp = -1*log10(pvalsC), adjustedP = p.adjust(pvalsC, method="fdr"))



#plasmid_assoc <- plasmid_assoc[c(1:11, 13:24, 26),]

# parse plasmid names for plotting / output
# generate plots of plasmid-specific ORs and a volcano plot

pOR <- plasmid_assoc %>% ggplot(aes(x = plasmid, y = OR)) + geom_point() + 
  geom_errorbar(aes(ymin = conflower, ymax = confupper)) + 
  coord_cartesian(ylim = c(0,5)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90)) + 
  labs(y = "Odds Ratio of Dissemination", x = "Plasmid Name")


pVolc <- plasmid_assoc %>% ggplot(aes(x = OR, y = logp, label = plasmid)) + 
  geom_point(alpha = 0.8) + 
  geom_text_repel(min.segment.length = 0) +
  theme_bw() + 
  labs(y = "-log10(p)", x = "Odds Ratio of Dissemination")
pVolc
plot_grid(pOR, pVolc, labels = c("C", "D"))
ggsave("../results/figures/paper/revision/Figure_4C-D.jpg", height = 4, width = 8)

plasmid_assoc$pval_adjusted_fdr <- p.adjust(plasmid_assoc$pval, method = "fdr")
write_tsv(plasmid_assoc, "../results/tables/submission/Supplemental_Table_4.tsv")


# Calculate plasmid-specific assocations w/dissemination for minimap analysis ####
## Conduct plasmid-specific presence/absence associations with dissemination
# this analysis uses the B31 minimap analysis
# create merged dataset
plas2 <- data.frame(plas2)
percent_occupancy = colSums(plas2) / nrow(plas2)
plas3 <- plas2[,percent_occupancy != 1]

plasmid_merge2 <- data.frame(plas3)
plasmid_merge2$label <- rownames(plas3)
plasmid_merge2 <- inner_join(plasmid_merge2, metadata_extended)

# conduct association tests
nplasmids <- ncol(plas3)
plasnames <- rep(NA, nplasmids)
pvals <- rep(NA, nplasmids)
estim <- rep(NA, nplasmids)
lowerint <- rep(NA, nplasmids)
upperint <- rep(NA, nplasmids)
for(i in c(1:nplasmids)){
  test_tmp <- fisher.test(plasmid_merge2[,i], plasmid_merge2$dissem_bin)
  plasnames[i] <- colnames(plasmid_merge2)[i]
  pvals[i] <- test_tmp$p.value
  estim[i] <- test_tmp$estimate
  lowerint[i] <- test_tmp$conf.int[1]
  upperint[i] <- test_tmp$conf.int[2]
  #print(paste(i, test_tmp$p.value, " ", colnames(plasmid_merge)[i]))   
}
plasmid_assoc2 <- data.frame(plasmid = plasnames, pval = pvals, OR = estim, conflower = lowerint, confupper = upperint, logp = -1*log10(pvals),
                             adjustedP = p.adjust(pvals, method="fdr"))

# parse plasmid names for plotting / output
#plasmid_assoc2$plasmid_name <- sapply(plasmid_assoc$plasmid, function(x) tail(strsplit(x, "_")[[1]], 1))
#plasmid_assoc2$plasmid_name <- gsub("ID.23.B31", "", plasmid_assoc$plasmid_name)

# generate plots of plasmid-specific ORs and a volcano plot

pOR2 <- plasmid_assoc2 %>% ggplot(aes(x = plasnames, y = OR)) + geom_point() + 
  geom_errorbar(aes(ymin = conflower, ymax = confupper)) + 
  coord_cartesian(ylim = c(0,5)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90)) + 
  labs(y = "Odds Ratio of Dissemination", x = "Plasmid Name")


pVolc2 <- plasmid_assoc2 %>% ggplot(aes(x = OR, y = logp, label = plasnames)) + 
  geom_point(alpha = 0.8) + 
  geom_text_repel(min.segment.length = 0) +
  theme_bw() + 
  labs(y = "-log10(p)", x = "Odds Ratio of Dissemination") + 
  coord_cartesian(xlim = c(0,3))

plot_grid(pOR2, pVolc2, labels = c("B", "C"))
ggsave("../results/figures/paper/revision/Figure_S5B_C.jpg", height = 4.5, width = 9)

plasmid_assoc2$pval_adjusted_fdr <- p.adjust(plasmid_assoc2$pval, method = "fdr")
write_tsv(plasmid_assoc2, "../results/tables/submission/Supplemental_Table_5.tsv")

# gene-specific associations ####
gt_merge <- data.frame(t(gt.annot))
gt_merge <- gt_merge[,colSums(gt_merge) < 297 & colSums(gt_merge) > 1]
ngenes <- ncol(gt_merge)
gt_merge$label <- rownames(gt_merge)
gt_merge <- inner_join(gt_merge, metadata_extended)

##
gene_p <- rep(NA, ngenes)
gene_eff <- rep(NA, ngenes)
l_int <- rep(NA, ngenes)
u_int <- rep(NA, ngenes)
for(i in 1:ngenes){
  test_stat <- fisher.test(gt_merge[,i], gt_merge$dissem_bin)
  gene_p[i] <- test_stat$p.value
  gene_eff[i] <- test_stat$estimate
  l_int[i] <- test_stat$conf.int[1]
  u_int[i] <- test_stat$conf.int[2]
}
##
gene_assoc <- data.frame(gene = colnames(gt_merge)[1:ngenes], pval = gene_p, OR = gene_eff, logp = -1*log10(gene_p), conf_lower = l_int, conf_upper = u_int,
                         adjustedP = p.adjust(gene_p, method="fdr")) %>%
  mutate(Annotation = str_split_fixed(gene, "\\.", 2)[,1], groupID = sapply(gene, function(x) tail(strsplit(x, "\\.")[[1]],1))) %>% mutate(ortho = paste(groupID, Annotation, sep="-")) %>% 
  mutate(desc = gsub("Borrelia\\.burgdorferi", "", gene)) %>% mutate(desc = gsub("\\.", " ", desc))

##
annot2 <- annot %>% rename("groupID" = Gene, "Annotation" = Locus) %>% mutate(ortho = paste(groupID, Annotation, sep="-"))

gene_assoc2 <- left_join(gene_assoc, annot2, by = c("ortho"))

gene_assoc2 %>% ggplot(aes(x = OR, y = logp, label = ortho, color = Localization)) + geom_text() + scale_x_log10()

# Figure 5C ####

pGAS <- gene_assoc2 %>% filter(Localization == "S", pval < 0.15) %>% ggplot(aes( x = desc, y = OR)) + geom_point() + geom_errorbar(aes(ymin = conf_lower, ymax = conf_upper)) + theme_bw() + #coord_cartesian(ylim = c(0,5)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, size = 7)) + scale_y_log10() + labs(x = "", y = "OR of Dissemination")
plot_grid(pGAS, labels = c("C"))
ggsave("../results/figures/paper/revision/Figure_5C.jpg", height = 6, width = 8)



heatmap(as.matrix(gt_merge[,1:ngenes]))

# Figure 7D ####
seer <- read_csv("../results/tables/submission/Supplemental_Table_6.csv")
ga3 <- gene_assoc2 %>% select(groupID.x, pval, OR) %>% rename("group" = groupID.x)
ga3 <- left_join(ga3, seer)

ga3 %>% ggplot(aes(x = pval, y = `filter-pvalue`, label = B31_Annotation)) + geom_point() + scale_x_log10() + scale_y_log10()

ga3 %>% #filter(af > 0.05 & af < 0.95) %>% 
  mutate(sublabel = ifelse(`lrt-pvalue` < 0.05, Locus, "")) %>%
  ggplot(aes(y = -1*log10(`lrt-pvalue`), x = exp(beta), label = sublabel)) + geom_text() + scale_x_log10()

pga3 <- ga3 %>% filter(`lrt-pvalue` < 0.1) %>% filter(af > 0.1 & af < 0.9) %>% 
  ggplot(aes( x = paste(Locus, Group, B31_Annotation, sep=" | "), y = exp(beta))) + geom_point() + geom_errorbar(aes(ymin = exp(beta - 1.96*`beta-std-err`), ymax = exp(beta + 1.96*`beta-std-err`))) + theme_bw() + coord_cartesian(ylim = c(0.01,10)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, size = 7)) + scale_y_log10() + labs(x = "", y = "OR of Dissemination")
plot_grid(pga3, labels = c("D"))
ggsave("../results/figures/paper/revision/Figure_7D.jpg", height = 6, width = 8)

# Figure S2A-F ####
## compare MCC and ML trees

beast_tree <- ggtree(b) %<+% metadata_extended + 
  #geom_nodelab(aes(label = round(posterior, 3), subset = posterior > 0.9), size = 3) + 
  geom_nodepoint(aes(color = posterior, subset = posterior > 0.9)) + 
  theme_tree2()

p1 <- ggtree(x) %<+% metadata_extended
p2 <- ggtree(b) %<+% metadata_extended

d1 <- p1$data
d2 <- p2$data

### reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree

p1$data$x <- p1$data$x
d2$x <- max(d2$x) - d2$x + max(d1$x) + 0.03
pp <- p1 + geom_tree(data=d2)
d1 <- p1$data

dd <- bind_rows(d1, d2) %>% filter(isTip == TRUE) %>% filter(!is.na(OspC_Type))
f2f3 <- pp + geom_line(aes(x, y, group=label, color=wgs_group), data=dd, alpha = 0.3) +  
  scale_color_discrete(name = "WGS Group")
US <- list(US = metadata_extended$label[metadata_extended$Region == "US Northeast" | metadata_extended$Region == "US Midwest"])
Midwest <- list(US = metadata_extended$label[metadata_extended$Region == "US Northeast" | metadata_extended$Region == "EU Slovenia"])
eu_tree <- groupOTU(b,  US)
mw_tree <- groupOTU(b, Midwest)

eu <- ggtree(eu_tree, aes(color = group)) + 
  scale_color_discrete(labels = c("EU", "US"), name = "") + 
  geom_nodelab(aes(label = round(posterior, 3), subset = posterior > 0.9), color = "black", size = 3) + 
  theme_tree2()

mw <- ggtree(mw_tree, aes(color = group)) + 
  geom_nodelab(aes(label = round(posterior, 3), subset = posterior > 0.9), color = "black" ,size = 3) + 
  scale_color_discrete(labels = c("Midwest", "Non-Midwest"), name = "") + 
  theme_tree2()

## compare MCC and ML trees


b2 <- read.beast("../results/BEAST/GTRG_constant_pop_gamma_dist_100M/mcc.tree")
metadata_extended2 <- metadata_extended %>%
  mutate(label = paste(label, Location, Year, sep= "|"))
beast_tree2 <- ggtree(b2) %<+% metadata_extended2 + 
  #geom_nodelab(aes(label = round(posterior, 3), subset = posterior > 0.9), size = 3) + 
  #geom_nodepoint(aes(color = posterior, subset = posterior > 0.9)) + 
  geom_tippoint(aes(color=Region, subset=!is.na(Region)), size = 1, alpha = 0.8) + 
  theme_tree2() + 
  geom_range(range='height_0.95_HPD', color='gray', alpha=.5, size=1) + 
  coord_cartesian(x = c(0, 3000000)) + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +
  labs(x = "Time (years)")

beast_trace <- readLog("../results/BEAST/GTRG_const_gamma_2e9_20M/core_gene_alignment_BEAST.aln.log", burnin = 0.3575)
bt <- data.frame(beast_trace)
bt$iteration <- rownames(bt)
bt <- bt[,c("iteration","treeModel.rootHeight", "tmrca.Slovenia.", "tmrca.Midwest.", "tmrca.Northeast.")]
names(bt) <- c("iteration", "All", "Slovenia", "Midwest", "Northeast")
bt_reshaped <- pivot_longer(bt, cols=-c(iteration), names_to = "Region", values_to = "TMRCA")

beast_tmrca <- ggplot(bt_reshaped, aes(x = Region, y = TMRCA, fill = Region)) + 
  geom_violin() + coord_cartesian(y = c(0, 5e6)) + 
  geom_boxplot(width = 0.1) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) + 
  labs(y = "TMRCA (Years)")

PG1 <- plot_grid(beast_tree, f2f3, eu, mw, beast_tree2, beast_tmrca, labels = c("A", "B", "C", "D", "E", "F"), rel_widths = c(1.75,2), nrow = 3)

# make plots of TMRCA under other models
btReshape <- function(BT, clockrate){
  BT <- data.frame(BT)
  BT$iteration <- rownames(BT)
  BT <- BT[,c("iteration","treeModel.rootHeight", "tmrca.Slovenia.", "tmrca.Midwest.", "tmrca.Northeast.")]
  names(BT) <- c("iteration", "All", "Slovenia", "Midwest", "Northeast")
  BT_r <- pivot_longer(BT, cols=-c(iteration), names_to = "Region", values_to = "TMRCA")
  BT_r$Clock_rate <- clockrate
  BT_r
}

bt8 <- readLog("../results/BEAST/GTRG_fixed_clock_1e8_20M/core_gene_alignment_BEAST.aln.log", burnin = 0.1)
BT8 <- btReshape(bt8, 1e-8)
bt9 <- readLog("../results/BEAST/GTRG_fixed_clock_1e9_20M/core_gene_alignment_BEAST.aln.log", burnin = 0.1)
BT9 <- btReshape(bt9, 1e-9)
bt10 <- readLog("../results/BEAST/GTRG_fixed_clock_1e10_20M/core_gene_alignment_BEAST.aln.log", burnin = 0.1)
BT10 <- btReshape(bt10, 1e-10)

BT_all1 <- rbind(BT8, BT9)
BT_all <- rbind(BT_all1, BT10)

pg2 <- ggplot(BT_all, aes(x = Region, y = TMRCA, fill = Region)) + 
  geom_violin() + #coord_cartesian(y = c(0, 5e6)) + 
  geom_boxplot(width = 0.1) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.position = "none") + 
  labs(y = "TMRCA (Years)") + 
  facet_wrap(~Clock_rate, scales = "free") 

PG2 <- plot_grid(pg2, labels = c("G"))
plot_grid(PG1, PG2, nrow = 2, rel_heights = c(3,1))
ggsave("../results/figures/paper/revision/Figure_S2A-G.jpg", height = 12, width = 10)


# Scrap ####