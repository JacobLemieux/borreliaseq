from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import argparse

parser = argparse.ArgumentParser(description = 'Blast pan-genome sequences from a multi fasta file')
parser.add_argument('-i', '--infile', help = 'input file', required=True)
parser.add_argument('-o', '--outfile', help = 'output filename', required=True)
parser.add_argument('-a', '--associationfile', help = 'file with p values for each ortholog from a panGWAS')
parser.add_argument('-b', '--blasttype', help = 'blastn or blastp', required=True)
parser.add_argument('-d', '--blastdb', help = 'blast DB, e.g. blast/B_burgdorfer_B31.fa, blast/PBr')
args = vars(parser.parse_args())

thresh = 0.000001

ortholog_dict = {}
f1 = open(args['associationfile'], "r")
for line in f1.readlines()[1:]:
	ortholog_dict[line.split()[0]] = float(line.split()[2])


blastdbpath = "blast/" + args['blastdb']

scaffold_dict = SeqIO.to_dict(SeqIO.parse(blastdbpath+".fna", "fasta"))
result_list = []

f = open(args['outfile']+"blast_annotations.out", "w")
for seq_record in SeqIO.parse(args['infile'], "fasta"):
	if ortholog_dict[seq_record.description.split(" ")[1]] < thresh:
		print("blasting", seq_record.id,"...")
		outfile_name = "tmp.xml"
		SeqIO.write(seq_record, "tmp.fasta", "fasta")
		#if args["blasttype"] == "blastp":
		#	blastp_cline = NcbiblastpCommandline(query="tmp.fasta", db="blast/pfam32_protein.fa", outfmt=5, out=outfile_name)
		if args["blasttype"] == "blastn":
			blastp_cline = NcbiblastnCommandline(query="tmp.fasta", db=blastdbpath, outfmt=5, out=outfile_name)
		stdout, stderr = blastp_cline()
		result_handle = open(outfile_name)
		blast_record = NCBIXML.read(result_handle)
		#print(blast_record)
		counter = 0
		new_scaffolds = []
		for alignment in blast_record.alignments:
			if counter < 1:
				print("Sequence name:", alignment.title)
				result_list.append(alignment.title.split(" ")[1])
				print(result_list[-1])
				out_line1 = "_".join(seq_record.id.split("_")[0:4]) # change to 2 if working with genomes from genbank
				out_line2 = "_".join(alignment.title.split("_")[2:4]).replace("No definition line", "")
				out_line_concat = out_line1 + "\t" + out_line2
				for hsp in alignment.hsps:
					print(hsp.positives,len(hsp.sbjct))
					out_line_concat = out_line_concat + "\t"  + scaffold_dict[result_list[-1]].description + "\t" + str(hsp.sbjct_start)
			#print("Sequence: ", alignment.hit)
				out_line_concat = out_line_concat + "\n"
				f.write(out_line_concat)
			counter = counter + 1
f.close()

# print only unique contigs
identifiers = set(result_list)
out_list = [scaffold_dict[x] for x in identifiers]
out_list.append(SeqIO.read("genome/Bb_chromosome_only.fna", "fasta"))
SeqIO.write(out_list, args['outfile'], "fasta")