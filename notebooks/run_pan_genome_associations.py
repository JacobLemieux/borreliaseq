import pandas as pd
import os
import numpy as np
import seaborn as sns
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import chi2_contingency

os.chdir(os.path.expanduser("~") + "/Dropbox/Tick/borrelia/borreliaseq/notebooks")
metadata = pd.read_csv("../results/LymeSeq_SampleTrack - Renaming Scheme Table for Jupyter.2023-02-14.csv", na_values="ND")
metadata["MEM_Binary"] = metadata['Multiple_EM(Y/N)'].replace(['N','Y'], [0,1])
metadata["Disseminated_Binary"] = metadata['Disseminated'].replace(['L','D'], [0,1])
metadata["OspC_A"] = [1 if i == "A"  else 0 for i in metadata["OspC_Type"]]
metadata["RST_1"] = [1 if i == 1  else 0 for i in metadata["RST_Type"]]
metadata["OspC_K"] = [1 if i == "K"  else 0 for i in metadata["OspC_Type"]]

pheno_output = metadata[["Rename_A", "Disseminated_Binary"]]
pheno_output.rename(columns = {"Rename_A": "Sample_name"}, inplace = True)
pheno_output["Sample_name"].replace('', np.nan, inplace = True)
pheno_output.dropna(subset = ["Sample_name"], inplace=True)
pheno_output.to_csv("../results/dissemination_phenotype.tsv",  float_format='%.0f', index = False, sep="\t")

pheno_output2 = metadata[["Rename_A", "OspC_A"]]
pheno_output2.rename(columns = {"Rename_A": "Sample_name"}, inplace = True)
pheno_output2["Sample_name"].replace('', np.nan, inplace = True)
pheno_output2.dropna(subset = ["Sample_name"], inplace=True)
pheno_output2.to_csv("../results/OspC_A_phenotype.tsv",  float_format='%.0f', index = False, sep="\t")

pheno_output3 = metadata[["Rename_A", "MEM_Binary"]]
pheno_output3.rename(columns = {"Rename_A": "Sample_name"}, inplace = True)
pheno_output3["Sample_name"].replace('', np.nan, inplace = True)
pheno_output3.dropna(subset = ["Sample_name"], inplace=True)
pheno_output3.to_csv("../results/MEM_phenotype.tsv",  float_format='%.0f', index = False, sep="\t")

pheno_output4 = metadata[["Rename_A", "RST_1"]]
pheno_output4.rename(columns = {"Rename_A": "Sample_name"}, inplace = True)
pheno_output4["Sample_name"].replace('', np.nan, inplace = True)
pheno_output4.dropna(subset = ["Sample_name"], inplace=True)
pheno_output4.to_csv("../results/RST_1_phenotype.tsv",  float_format='%.0f', index = False, sep="\t")

pheno_output5 = metadata[["Rename_A", "OspC_K"]]
pheno_output5.rename(columns = {"Rename_A": "Sample_name"}, inplace = True)
pheno_output5["Sample_name"].replace('', np.nan, inplace = True)
pheno_output5.dropna(subset = ["Sample_name"], inplace=True)
pheno_output5.to_csv("../results/OspC_K_phenotype.tsv",  float_format='%.0f', index = False, sep="\t")


## create distance matrix from tree
#os.system("python ../scripts/phylogeny_distances.py ../results/roary-results/core_gene_alignment_midpoint_rooted.nwk > ../results/roary-results/distance_matrix_from_tree.tsv")

## Create similarity matrix from tree
#os.system("python ../scripts/phylogeny_distances.py --lmm ../results/roary-results/core_gene_alignment_midpoint_rooted.nwk > ../results/roary-results/similarity_matrix_from_tree.tsv")

## run the analysis without lineage correction
#os.system("pyseer --pres ../results/roary-results/gene_presence_absence.Rtab --phenotypes ../results/dissemination_phenotype.tsv --no-distances > ../results/pyseer_no_distances_dissemination_results.txt")

## run pyseer with default lineage correction
#os.system("pyseer --pres ../results/roary-results/gene_presence_absence.Rtab --phenotypes ../results/dissemination_phenotype.tsv --distances ../results/roary-results/distance_matrix_from_tree.tsv > ../results/dissemination_results_pyseer_default.txt")

## run pyseer with default lineage correction and max dimension
#os.system("pyseer --pres ../results/roary-results/gene_presence_absence.Rtab --phenotypes ../results/dissemination_phenotype.tsv --distances ../results/roary-results/distance_matrix_from_tree.tsv --max-dimension 6 > ../results/dissemination_results_pyseer_default_max_dim_6.txt")

# run pyseer with LMM and similarity matrix
#os.system("pyseer --pres ../results/roary-results/gene_presence_absence.Rtab --phenotypes ../results/dissemination_phenotype.tsv --lmm --similarity ../results/roary-results/similarity_matrix_from_tree.tsv > ../results/dissemination_results_lmm.txt")

### This is the model reported in the paper
if False:
  # run with lineage effects
  os.system("pyseer --pres ../results/roary-results/gene_presence_absence.Rtab --phenotypes ../results/dissemination_phenotype.tsv --lineage --lineage-file ../results/roary-results/inferred_lineages.txt --distances ../results/roary-results/distance_matrix_from_tree.tsv > ../results/dissemination_results_lineage.txt")
  
  ## run the analysis without lineage correction to identify OspC Type A-associated cases
  os.system("pyseer --pres ../results/roary-results/gene_presence_absence.Rtab --phenotypes ../results/OspC_A_phenotype.tsv --no-distances > ../results/pyseer_no_distances_OspC_A_results.txt")
  
  ## run the analysis without lineage correction to identify RST 1-associated cases
  os.system("pyseer --pres ../results/roary-results/gene_presence_absence.Rtab --phenotypes ../results/RST_1_phenotype.tsv --no-distances > ../results/pyseer_no_distances_RST_1_results.txt")
  
  ## run the analysis without lineage correction to identify OspC Type K-associated cases
  os.system("pyseer --pres ../results/roary-results/gene_presence_absence.Rtab --phenotypes ../results/OspC_K_phenotype.tsv --no-distances > ../results/pyseer_no_distances_OspC_K_results.txt")
  ## run the analysis with lineage correction to identify MEM-associated cases
  os.system("pyseer --pres ../results/roary-results/gene_presence_absence.Rtab --phenotypes ../results/MEM_phenotype.tsv --lineage --lineage-file ../results/roary-results/inferred_lineages.txt --distances ../results/roary-results/distance_matrix_from_tree.tsv > ../results/pyseer_no_distances_MEM_results.txt")

if False:
  ## Run Scoary (Requires phenotypes as integer rather than float)
  os.system("/opt/anaconda3/bin/scoary -g ../results/roary-results/gene_presence_absence.csv -t ../results/dissemination_phenotype_scoary.csv -o ../results/scoary/")

# concatenate B burgdorferi genome, map pan genome elements to reference using minimap2.
if False:
    seqsum = SeqRecord(Seq(""))
    for seqrecord in SeqIO.parse("../blast/B_burgdorferi_B31.fa.fna", "fasta"):
        seqsum.seq = seqsum.seq + seqrecord.seq
    seqsum.id = "Concatenated B31"
    SeqIO.write(seqsum, "../genome/B_burgdorferi_B31_concatenated.fa.fna", format = "fasta")
    os.system("minimap2 -x asm5 ../genome/B_burgdorferi_B31_concatenated.fa.fna ../results/roary-results/pan_genome_reference.fa > ../alignment/pan_genome_alignment.paf")
if False:
    os.system("minimap2 -x asm5 ../blast/B_burgdorferi_B31.fa.fna ../results/roary-results/pan_genome_reference.fa > ../alignment/pan_genome_alignment_by_chromosome.paf")


# obtain chromosome lengths for B burgdorferi reference genome
seqlengths = [len(seqrecord)for seqrecord in SeqIO.parse("../blast/B_burgdorferi_B31.fa.fna", "fasta")]
seqnames = [seqrecord.id for seqrecord in SeqIO.parse("../blast/B_burgdorferi_B31.fa.fna", "fasta")]
chr_boundary = np.cumsum(seqlengths)[-2]
chr_coords = pd.DataFrame({'Replicon': seqnames, 
                           'Replicon_Offset': np.pad(np.cumsum(seqlengths)[:-1], (1,0)),
                           'Replicon_End': np.cumsum(seqlengths)
                          } )

# Read in minimap2 alignment results
alignment_results = pd.read_csv("../alignment/pan_genome_alignment_by_chromosome.paf", sep="\t", header=None)
alignment_results.columns = ["ID", "qlength", "qstart", "qend", "strand", "refname", "rlength","rstart", "rend", "rmatch", "rexact", "qual", "a", "b", "c", "d", "e", "f"]
for index, row in alignment_results.iterrows():
    alignment_results.loc[index,'Replicon_Coordinate'] = int(alignment_results["rstart"][index] + chr_coords.loc[chr_coords["Replicon"] == row["refname"]]["Replicon_Offset"])
alignment_results = alignment_results.rename(columns={"refname": "Replicon"})

# obtain mapping between sequence IDs and group IDs and add to alignment results
pan_genome_seqs = {}
for seqrecord in SeqIO.parse("../results/roary-results/pan_genome_reference.fa", "fasta"):
    pan_genome_seqs[seqrecord.id] = seqrecord.description
pan_genome_IDs = pd.DataFrame.from_dict(pan_genome_seqs, orient='index', columns=["ID"])
pan_genome_IDs[["ID", "group"]] = pan_genome_IDs["ID"].str.split(" ", 1, expand=True)


# merge pan_genome IDs
alignment_results = pan_genome_IDs.merge(alignment_results, on = "ID")

# Read in gene annotation from Roary
annotation = pd.read_csv("../results/roary-results/gene_presence_absence.csv")[["Gene", "Annotation"]]
annotation = annotation.rename(columns = {"Gene":"group"})
alignment_results = alignment_results.merge(annotation)
more_annotation = pd.read_csv("../results/roary-results/pan_genome_with_annotation.csv")
alignment_results = alignment_results.merge(more_annotation)
lipo_annotation = pd.read_csv("../results/annotation/Dowdell_lipoproteins.csv")
alignment_results = alignment_results.merge(lipo_annotation, how = "left")
alignment_results["Localization2"] = alignment_results["Localization"].replace(np.nan, "unknown", regex=True)
alignment_results.to_csv("../results/annotation/compiled_annotation.csv")

# Read in results of pyseer
unadj = pd.read_csv("../results/pyseer_no_distances_dissemination_results.txt", sep="\t", index_col=False)
unadj = unadj.rename(columns={"variant": "group"})
unadj = unadj.merge(alignment_results, on = "group")
unadj["-log10 P"] = -1 * np.log10(unadj["filter-pvalue"])

## Pyseer with default lineage correction, no max dimension. This is deprecated.
#lin_adj = pd.read_csv("../results/dissemination_results_pyseer_default.txt", sep="\t", index_col=False)
#lin_adj = lin_adj.rename(columns={"variant": "group"})
#lin_adj["-log10 P"] = -1 * np.log10(lin_adj["lrt-pvalue"])
#lin_adj = lin_adj.merge(alignment_results, on = "group")

## Pyseer with lineage correction, max dimension 6. This is deprecated.
#lin_adj_md6 = pd.read_csv("../results/dissemination_results_pyseer_default_max_dim_6.txt", sep="\t", index_col=False)
#lin_adj_md6 = lin_adj_md6.rename(columns={"variant": "group"})
#lin_adj_md6["-log10 P"] = -1 * np.log10(lin_adj_md6["lrt-pvalue"])
#lin_adj_md6 = lin_adj_md6.merge(alignment_results, on = "group")

## Pyseer with lmm option (bugWAS model). Very similar to lineage correction, and not reported in the paper.
#lmm = pd.read_csv("../results/dissemination_results_lmm.txt", sep="\t", index_col=False)
#lmm = lmm.rename(columns={"variant": "group"})
#lmm["-log10 P"] = -1 * np.log10(lmm["lrt-pvalue"])
#lmm = lmm.merge(alignment_results, on = "group")

scoary = pd.read_csv("../results/scoary/_20_05_2023_1725.results.csv")
scoary = scoary.rename(columns={"Gene":"group", "Annotation":"Scoary_Annotation"})
scoary = scoary.merge(alignment_results, on = "group")

# OspC type A associated. Reported in paper.
ospA = pd.read_csv("../results/pyseer_no_distances_OspC_A_results.txt", sep="\t", index_col=False)
ospA = ospA.rename(columns={"variant": "group"})
ospA["-log10 P"] = -1 * np.log10(ospA["lrt-pvalue"])
ospA2 = ospA.merge(more_annotation, left_on = "group", right_on = "Group")
ospA = ospA.merge(alignment_results, on = "group")

# OspC type K associated. Reported in paper.
ospK = pd.read_csv("../results/pyseer_no_distances_OspC_K_results.txt", sep="\t", index_col=False)
ospK = ospK.rename(columns={"variant": "group"})
ospK["-log10 P"] = -1 * np.log10(ospK["lrt-pvalue"])
ospK2 = ospK.merge(more_annotation, left_on = "group", right_on = "Group")
ospK = ospK.merge(alignment_results, on = "group")

# OspC type K associated. Reported in paper.
RST = pd.read_csv("../results/pyseer_no_distances_RST_1_results.txt", sep="\t", index_col=False)
RST = RST.rename(columns={"variant": "group"})
RST["-log10 P"] = -1 * np.log10(RST["lrt-pvalue"])
RST2 = RST.merge(more_annotation, left_on = "group", right_on = "Group")
RST = RST.merge(alignment_results, on = "group")


# Pyseer with lineage correction and output of inferred lineages to results/roary_results/inferred_lineages.txt. Reported in paper.
lineage = pd.read_csv("../results/dissemination_results_lineage.txt", sep="\t", index_col=False)
lineage = lineage.rename(columns={"variant": "group"})
lineage["-log10 P"] = -1 * np.log10(lineage["lrt-pvalue"])
lineage2 = lineage.merge(more_annotation, left_on = "group", right_on = "Group")
lineage = lineage.merge(alignment_results, on = "group")

lineage2 = lineage2.merge(lipo_annotation, on = "Locus", how = 'left')

fig, axes = plt.subplots(3,1)
ospA_lipo = ospA[ospA["Localization"].isin(["S", "P-IM", "P-OM"])]
sns.scatterplot(ax = axes[0],x = ospA_lipo["Replicon_Coordinate"], y = -1*np.log10(ospA_lipo["filter-pvalue"]), hue = ospA_lipo["Localization"])
plt.legend(loc='upper right')
plt.xlabel("Genomic Coordinate")
plt.ylabel("-log10(P)")

ospK_lipo = ospK[ospK["Localization"].isin(["S", "P-IM", "P-OM"])]
sns.scatterplot(ax = axes[1],x = ospK_lipo["Replicon_Coordinate"], y = -1*np.log10(ospK_lipo["filter-pvalue"]), hue = ospK_lipo["Localization"])
plt.legend(loc='upper right')
plt.xlabel("Genomic Coordinate")
plt.ylabel("-log10(P)")

RST_lipo = RST[RST["Localization"].isin(["S", "P-IM", "P-OM"])]
sns.scatterplot(ax = axes[2],x = RST_lipo["Replicon_Coordinate"], y = -1*np.log10(RST_lipo["filter-pvalue"]), hue = RST_lipo["Localization"])
plt.legend(loc='upper right')
plt.xlabel("Genomic Coordinate")
plt.ylabel("-log10(P)")

plt.savefig("../results/figures/paper/Figure_S7.png")
plt.close()

#ospA["Linked"] = [True if i < 1e-3 else False for i in ospA["filter-pvalue"]]
#pd.crosstab(ospA["Localization2"], ospA["Linked"])
xvals = [0, 200000, 400000, 600000, 800000, 1000000, 1200000, 1400000]
xlabs = ['0', '200', '400', '600', '800', '1,000', '1,200', '1,400']

plt.figure(figsize=(14,2.5))
p6_1 = sns.scatterplot(x = lineage["Replicon_Coordinate"], y = -1*np.log10(lineage["filter-pvalue"]), hue = lineage["Replicon"])
#for location in np.cumsum(seqlengths):
    #p6_1.axvline(location, alpha = 0.5)
plt.xticks(xvals,xlabs, fontdict={'fontsize':16})
plt.xlabel("Concatenated Genomic Coordinate (Kilobases)", fontdict={'fontsize':16})
plt.ylabel("-log10(P)")
plt.title("ORF-Associations with Dissemination, no Lineage Correction", fontdict={'fontsize':18})
plt.legend(bbox_to_anchor=(1.2, 1.1))
p6_1.legend_.remove()
plt.savefig("../results/figures/paper/revision/Figure_7A.pdf", bbox_inches="tight")
plt.close()

plt.figure(figsize=(14,8.5))
p6_1 = sns.scatterplot(x = lineage["Replicon_Coordinate"], y = -1*np.log10(lineage["filter-pvalue"]), hue = lineage["Replicon"])
#for location in np.cumsum(seqlengths):
    #p6_1.axvline(location, alpha = 0.5)
plt.xticks(xvals,xlabs, fontdict={'fontsize':16})
plt.xlabel("Concatenated Genomic Coordinate (Kilobases)", fontdict={'fontsize':16})
plt.ylabel("-log10(P)")
plt.legend(bbox_to_anchor=(1.12, 0.9))
p6_2.legend_.remove()
plt.savefig("../results/figures/paper/Figure_7A_with_legend.pdf", bbox_inches="tight")
plt.close()

plt.figure(figsize=(14,2.5))
p6_2 = sns.scatterplot(x = lineage["Replicon_Coordinate"], y = -1*np.log10(lineage["lrt-pvalue"]), hue = lineage["Replicon"])
plt.xticks(xvals,xlabs, fontdict={'fontsize':16})
plt.xlabel("Concatenated Genomic Coordinate (Kilobases)", fontdict={'fontsize':16})
plt.ylabel("-log10(P)")
plt.legend(bbox_to_anchor=(1.15, 0.8))
p6_2.legend_.remove()
plt.title("ORF-Associations with Dissemination, with Lineage Correction", fontdict={'fontsize':18})
plt.savefig("../results/figures/paper/revision/Figure_7B.pdf", bbox_inches="tight")
plt.close()

plt.figure(figsize=(14,2.5))
p6_2alt = sns.scatterplot(x = scoary["Replicon_Coordinate"], y = -1*np.log10(scoary["Naive_p"]), hue = scoary["Replicon"])
plt.xticks(xvals,xlabs, fontdict={'fontsize':16})
plt.xlabel("Concatenated Genomic Coordinate (Kilobases)", fontdict={'fontsize':16})
plt.ylabel("-log10(P)")
plt.legend(bbox_to_anchor=(1.15, 0.8))
#p6_2alt.legend_.remove()
plt.title("ORF-Associations with Dissemination, with Lineage Correction", fontdict={'fontsize':18})
plt.savefig("../results/figures/paper/revision/Figure_7B_alt.pdf", bbox_inches="tight")
plt.close()



plt.figure(figsize=(14,2.5))
lin_sub = lineage[lineage["lineage"].isin(["MDS8", "MDS10"])]
p6_21 = sns.scatterplot(x = lin_sub["Replicon_Coordinate"], y = -1*np.log10(lin_sub["filter-pvalue"]), hue = lin_sub["lineage"])
plt.xticks(xvals,xlabs, fontdict={'fontsize':16})
plt.xlabel("Concatenated Genomic Coordinate (Kilobases)", fontdict={'fontsize':16})
p6_21.legend_.remove()
plt.ylabel("-log10(P)")
plt.legend(bbox_to_anchor=(0.9, 0.9))
plt.title("ORF-Associations with Dissemination, by Lineage", fontdict={'fontsize':18})
plt.savefig("../results/figures/paper/revision/Figure_7C.pdf", bbox_inches="tight")
plt.close()

plt.figure(figsize=(14,2.5))
p6_3 = sns.scatterplot(x = ospA["Replicon_Coordinate"], y = -1*np.log10(ospA["lrt-pvalue"]), hue = ospA["Replicon"])
plt.xticks(xvals,xlabs, fontdict={'fontsize':16})
plt.xlabel("Concatenated Genomic Coordinate (Kilobases)", fontdict={'fontsize':16})
plt.ylabel("-log10(P)")
plt.legend(bbox_to_anchor=(1.15, 0.8))
p6_3.legend_.remove()
plt.title("ORF-Associations with OspC Type A Genotype", fontdict={'fontsize':18})
plt.savefig("../results/figures/paper/revision/Figure_8A.pdf", bbox_inches="tight")
plt.close()

plt.figure(figsize=(14,2.5))
p6_4 = sns.scatterplot(x = ospK["Replicon_Coordinate"], y = -1*np.log10(ospK["lrt-pvalue"]), hue = ospK["Replicon"])
plt.xticks(xvals,xlabs, fontdict={'fontsize':16})
plt.xlabel("Concatenated Genomic Coordinate (Kilobases)", fontdict={'fontsize':16})
plt.ylabel("-log10(P)")
plt.legend(bbox_to_anchor=(1.15, 0.8))
p6_4.legend_.remove()
plt.title("ORF-Associations with OspC Type K Genotype", fontdict={'fontsize':18})
plt.savefig("../results/figures/paper/revision/Figure_8B.pdf", bbox_inches="tight")
plt.close()

plt.figure(figsize=(14,2.5))
p6_5 = sns.scatterplot(x = RST["Replicon_Coordinate"], y = -1*np.log10(RST["lrt-pvalue"]), hue = RST["Replicon"])
plt.xticks(xvals,xlabs, fontdict={'fontsize':16})
plt.xlabel("Concatenated Genomic Coordinate (Kilobases)", fontdict={'fontsize':16})
plt.ylabel("-log10(P)")
plt.legend(bbox_to_anchor=(1.15, 0.8))
p6_5.legend_.remove()
plt.title("ORF-Associations with RST1 Genotype", fontdict={'fontsize':18})
plt.savefig("../results/figures/paper/revision/Figure_8C.pdf", bbox_inches="tight")
plt.close()



#sns.set(rc = {'figure.figsize':(14,2.5), 
#             'axes.facecolor':'white', 'figure.facecolor':'white'}
#       )

#lin_sub = lineage[lineage["Annotation"].str.contains("surface")]

## Output supplemental data files

# add multiple p correction using FDR
lineage3 = lineage2
filter_rejected, filter_corrected_p = sm.stats.fdrcorrection(lineage3["filter-pvalue"])
lrt_rejected, lrt_corrected_p = sm.stats.fdrcorrection(lineage3["lrt-pvalue"])
lineage3["filter-pvalue-fdr"] = filter_corrected_p
lineage3["lrt-pvalue-fdr"] = lrt_corrected_p

# supplemental data file 1
lineage3.sort_values(by = "lrt-pvalue").to_csv("../results/tables/submission/supplemental_data_file_1.csv")

# add multiple p correction using FDR
lineage4 = lineage2[lineage2["Localization"] == "S"]
filter_rejected, filter_corrected_p = sm.stats.fdrcorrection(lineage4["filter-pvalue"])
lrt_rejected, lrt_corrected_p = sm.stats.fdrcorrection(lineage4["lrt-pvalue"])
lineage4["filter-pvalue-fdr"] = filter_corrected_p
lineage4["lrt-pvalue-fdr"] = lrt_corrected_p

# supplemental data file 2
lineage4.sort_values(by = "lrt-pvalue").to_csv("../results/tables/submission/supplemental_data_file_2.csv")

# output table of genes associated with OspC Type A strains
ospA3 = ospA2
filter_rejected, filter_corrected_p = sm.stats.fdrcorrection(ospA3["filter-pvalue"])
lrt_rejected, lrt_corrected_p = sm.stats.fdrcorrection(ospA3["lrt-pvalue"])
ospA3["filter-pvalue-fdr"] = filter_corrected_p
ospA3["lrt-pvalue-fdr"] = lrt_corrected_p
ospA3.sort_values(by = "lrt-pvalue").to_csv("../results/tables/submission/supplemental_data_file_1b.csv")

# output table of genes associated with OspC Type A strains
ospK3 = ospK2
filter_rejected, filter_corrected_p = sm.stats.fdrcorrection(ospK3["filter-pvalue"])
lrt_rejected, lrt_corrected_p = sm.stats.fdrcorrection(ospK3["lrt-pvalue"])
ospK3["filter-pvalue-fdr"] = filter_corrected_p
ospK3["lrt-pvalue-fdr"] = lrt_corrected_p
ospK3.sort_values(by = "lrt-pvalue").to_csv("../results/tables/submission/supplemental_data_file_1c.csv")

# output table of genes associated with OspC Type A strains
RST3 = RST2
filter_rejected, filter_corrected_p = sm.stats.fdrcorrection(RST3["filter-pvalue"])
lrt_rejected, lrt_corrected_p = sm.stats.fdrcorrection(RST3["lrt-pvalue"])
RST3["filter-pvalue-fdr"] = filter_corrected_p
RST3["lrt-pvalue-fdr"] = lrt_corrected_p
RST3.sort_values(by = "lrt-pvalue").to_csv("../results/tables/submission/supplemental_data_file_1d.csv")
