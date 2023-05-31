from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import argparse

parser = argparse.ArgumentParser(description = 'Blast pfam32 sequences from a multi fasta file')
parser.add_argument('-i', '--infile', help = 'input file', required=True)
parser.add_argument('-o', '--outfile', help = 'output filename', required=True)
parser.add_argument('-b', '--blasttype', help = 'blastn or blastp or blastcasjens or blastcasjensnolp2811', required=True)
parser.add_argument('-g', '--generate_pseudoassemb', help = 'generate a pseudoassembly by adding -g yes')
parser.add_argument('-t', '--translate', help = 'translate input, no or yes', required=True)
args = vars(parser.parse_args())


#scaffold_dict = SeqIO.to_dict(SeqIO.parse("assembly/genbank/fasta/complete_plasmids.fna", "fasta"))
result_list = []

f = open(args['outfile']+"blast_annotations.out", "w")
for seq_record in SeqIO.parse(args['infile'], "fasta"):
	#print("blasting", seq_record.id,"...")
	if args['translate'] == "yes": 
		new_seq = SeqRecord(seq_record.seq.translate(), id = seq_record.id)
		#print(new_seq.seq)
	if args['translate'] == "no":
		new_seq = SeqRecord(seq_record.seq, id = seq_record.id)
		#print(new_seq.seq)
	outfile_name = "tmp.xml"
	SeqIO.write(new_seq, "tmp.fasta", "fasta")
	if args["blasttype"] == "blastp":
		blastp_cline = NcbiblastpCommandline(query="tmp.fasta", db="blast/pfam32_protein.fa", outfmt=5, out=outfile_name)
	if args["blasttype"] == "blastn":
		blastp_cline = NcbiblastnCommandline(query="tmp.fasta", db="blast/complete_plasmids.fna", outfmt=5, out=outfile_name)
	if args["blasttype"] == "blastcasjens":
		blastp_cline = NcbiblastpCommandline(query="tmp.fasta", db="blast/Casjens_pfam32_protein.fa", outfmt=5, out=outfile_name)
	if args["blasttype"] == "blastcasjensnolp2811":
		blastp_cline = NcbiblastpCommandline(query="tmp.fasta", db="blast/Casjens_pfam32_subset.fa", outfmt=5, out=outfile_name)
	stdout, stderr = blastp_cline()
	result_handle = open(outfile_name)
	blast_record = NCBIXML.read(result_handle)
	#print(blast_record)
	counter = 0
	new_scaffolds = []
	for alignment in blast_record.alignments:
		if counter < 3:
			print(alignment.hit_def + "\t" + alignment.hit_id + "\t" + alignment.accession + "\t" + seq_record.id)
			#print(dir(alignment))            
		if counter < 1:
			#print("Sequence name:", alignment.title)
			result_list.append(alignment.title.split(" ")[1])
			out_line1 = "_".join(seq_record.id.split("_")[0:4]) # change to 2 if working with genomes from genbank
			out_line2 = "_".join(alignment.title.split("_")[2:4]).replace("No definition line", "")
			out_line_concat = out_line1 + "\t" + out_line2 + "\n"
			#print("Sequence: ", alignment.hit)
			f.write(out_line_concat)
		counter = counter + 1
f.close()

# print only unique contigs
#identifiers = set(result_list)
#out_list = [scaffold_dict[x] for x in identifiers]

# append chromosome sequences (only do for pseudoassemblies)
if args['generate_pseudoassemb'] == "yes":
    out_list.append(SeqIO.read("genome/Bb_chromosome_only.fna", "fasta"))
print(result_list)
#SeqIO.write(result_list, args['outfile'], "fasta")