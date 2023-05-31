# Form pseudoassemblies from pfam32 sequences

import pandas as pd
metadata = pd.read_csv("results/LymeSeq_SampleTrack - Renaming Scheme Table for Jupyter.2023-01-14.csv")
metadata = metadata.dropna(subset=["Rename_A"])
#exclude_list = [] # ["UMA6", "ESI25", "ESI27", "URI45", "UCT95", "UCT124", "URI125", "URI126", "UMA128", "UNY152", "ESI288b", "ESI293", "ESI296", "ESI302", "ESI310", "ESI321"]
#SAMPLES = [item for item in metadata["Rename_A"].tolist() if item not in exclude_list]
SAMPLES = [item for item in metadata["Rename_A"].tolist()]
SAMPLES.append("B31ref")
THREADS = 10

rule all:
    input: 
        expand("results/annotation/short_read2/pfam32/{id}_pfam_subset.out", id=SAMPLES) 
        #expand("results/annotation/short_read2/pfam32/{id}_nucleotide.pfam32.outcontig_scaffold.fasta", id=SAMPLES)#, 
        #expand("results/ragoo/{id}/ragoo_output/{id}_ragoo.fasta", id=SAMPLES)

rule annotate_with_prokka:
    input:
        'results/assemblies/{id}_200.fa'
    output:
        'results/prokka/{id}/prokka_annotated.faa'
    threads: 10
    shell:
        "prokka {input} --proteins genome/GCF_000008685.2_ASM868v2_protein.faa --outdir results/prokka/{wildcards.id} --prefix prokka_annotated --cpus {threads} --rawproduct --force --genus Borrelia --kingdom Bacteria --centre X --compliant"

rule apply_pf32_hmm:
    input:
        "results/prokka/{id}/prokka_annotated.faa"
    output:
        "results/annotation/short_read2/{id}_pfam.out"
    shell:
        "hmmsearch --tblout {output} blast/pfam_full.hmm {input}"

rule extract_pf32_hits:
    input:
        "results/annotation/short_read2/{id}_pfam.out"
    output:
        "results/annotation/short_read2/pfam32/{id}_pfam_subset.out"
    shell:
        "python scripts/parse_hmmer_output.py -i results/prokka/{wildcards.id}/prokka_annotated.faa -m {input} -o {output}"


#rule blast_pfam32_sequences_from_short_read_assemblies:
#    input:
#        "results/annotation/short_read2/{id}_nucleotide.fa"
#    output:
#        "results/annotation/short_read2/pfam32/{id}_nucleotide.pfam32.outcontig_scaffold.fasta"
#    shell:
#        "python scripts/blast_pfam32_sequences.py -i {input} -o {output} -b blastn -t yes"

#rule execute_RaGOO: 
#	input:
#		"annotation/short_read/pfam32/{id}_nucleotide.pfam32.outcontig_scaffold.fasta"
#	output:
#		"ragoo/{id}/ragoo_output/{id}_ragoo.fasta"
#	run:
#		shell("~/bin/RaGOO-master/ragoo.py assembly/fasta/{wildcards.id}_assembly_spades.fa {input}")
#		shell("mv ragoo_output ragoo/{wildcards.id}/")
#		shell("mv ragoo/{wildcards.id}/ragoo_output/ragoo.fasta {output}")