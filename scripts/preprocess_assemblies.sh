## trimmomatic trimming ##
TrimmomaticPE -threads $(nproc) -phred33 id_R1_001.fastq.gz id_R2_001.fastq.gz id_R1.fastq.gz /dev/null id_R2.fastq.gz /dev/null ILLUMINACLIP:/usr/share/trimmomatic/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## spades assembly ##
python3.8 spades.py -k 21,33,55,77,99 --careful -1 id_R1.fastq.gz -2 id_R2.fastq.gz -o id_spades

## seqtk cutoff ##
seqtk seq -L 200 id_contigs.fa > id_contigs_200.fa  ## id_contigs.fa was generated from spades assembly ##

## assmebly review ##
grep -v ">" id_contigs_200.fa | wc | awk '{print $3-$1}'  ## assembly size ##
grep -c '^>' id_contigs_200.fa  ## contig numbers ##

## kraken2 classification, if needed ##
kraken2 --db Mini_Kraken2_10182017 --threads $(nproc) id_contigs_200.fa > id.kraken2	
cut -f2,3 id.kraken2 > id.krona.in;
ktImportTaxonomy id.krona.in -o id.krona.html

## assembly cleaning, if needed ##
mv id_contigs_200.fa id_contigs_ori.fa
cut -c 1- clean.list | xargs -n 1 samtools faidx id_contigs_ori.fa > id_contigs_200.fa ## clean.list was an ID list of contigs classified to Bb by kraken2 ##

## prokka annotation ##
prokka --outdir prokka/id --cpus $(nproc) --evalue 1e-10 --gram neg --kingdom Bacteria --genus Borrelia --species burgdorferi --prefix id --usegenus --addgenes --rfam id_contig_200.fa;

## roary core/pan-genome analysis ##
roary -f roary_results -e -n -v $(cat gff_files)  ## gff_files was a list of GFF files generated from prokka annotation ##

## kwip relatedness ##
khmer load-into-counting.py -N 1 -x 1e9 -k 31 -b -T $(nproc) -f -s tsv hash/id.ct.gz id_contigs_200.fa
kwip -t $(nproc) -k Bb.kern -d Bb.dist hash/*.ct.gz

## plots in RStudio ##
library(ctc)
library(fields)
dist = as.matrix(read.delim("Bb.dist")[, -1])
clust = hclust(as.dist(dist))
plot(clust, cex=0.55, main="Borrelia burgdorferi 299 isolates")
newick <- hc2Newick(clust)
write(newick, file="Bb.nwk")
fit <- cmdscale(dist, eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type='p')
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggstar)
library(treeio)
library(ggnewscale)
metadata <- read.csv("metadata.csv")
treeBb <- read.newick("Bb.nwk")
p <- ggtree(treeBb, layout = "fan", branch.length = "none", open.angle=10) %<+% metadata + geom_tiplab(aes(color = Location), size = 2.5, offset = 0.2) + scale_color_manual(values = c("#F8766D", "#C49A00", "#53B400", "#00B6EB", "#A58AFF", "#FB61D7"))
p
p2 <- p + new_scale_color() + new_scale_fill() + geom_tiplab2(aes(color = OspC, label = OspC), align = T, linetype = NA, size = 2.5, offset = 3.5, hjust = 0.5)
p2
p3 <- p2 + new_scale_color() + new_scale_fill() + geom_tiplab2(aes(color = RST, label = RST), align = T, linetype = NA, size = 2.5, offset = 4.5, hjust = 0.5)
p3
p4 <- p3 + new_scale_color() + new_scale_fill() + geom_tiplab2(aes(color = MLST, label = MLST), align = T, linetype = NA, size = 2.5, offset = 6, hjust = 0.5)
p4
