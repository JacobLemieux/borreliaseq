# Filter_assemblies using blast and coverage
# January 12 2020
# Jacob E. Lemieux

# pseudocode:
# read in assemblies by sequence
# blast each sequence against a custom blast db that consists of all Borrelia proteins
# parse coverage number
# hold sequences in a list, and output to new file those that have blast score exceeding certain level and minimum coverage

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import argparse
import matplotlib.pyplot as plt
import math
import numpy as np

# will have to loop through file, extract properties of the distribution of coverage and then loop through again to blast. 
# this will also speed up the blast. 

parser = argparse.ArgumentParser(description = 'Filter assemblies from a multi fasta file')
parser.add_argument('-i', '--infile', help = 'input file', required=True)
parser.add_argument('-o', '--outfile', help = 'output filename', required=True)
args = vars(parser.parse_args())
coverage_list = []
for seq_record in SeqIO.parse(args['infile'], "fasta"):
	#print("blasting", seq_record.id,"...")
	coverage = seq_record.id.split("_")[5]
	#print("Coverage",coverage, " ; log coverage", math.log(float(coverage)+0.0001))
	coverage_list.append(math.log(float(coverage)+0.0001))
	
mediancov = np.median(coverage_list)
print("Sample: ", len(coverage_list), " records;", args['infile'],"; median ln(coverage):", mediancov, "ln(coverage) standard deviation:", np.std(coverage_list))
#plt.hist(coverage_list)
#plt.show()
filtered_list = []
for seq_record in SeqIO.parse(args['infile'], "fasta"):
	#print("processing", seq_record.id,"...")
	coverage = seq_record.id.split("_")[5]
	if (mediancov - math.log(float(coverage)+0.0001)) < 2.3: # remove things with 10-fold less coverage than median
		#print("retaining", seq_record.id)
		outfile_name = args['infile']+".tmp.xml"
		SeqIO.write(seq_record, args['infile']+".tmp.fasta", "fasta")
		blastn_cline = NcbiblastnCommandline(query=args['infile']+".tmp.fasta", db="blast/combined.fna", outfmt=5, out=outfile_name)
		stdout, stderr = blastn_cline()
		result_handle = open(outfile_name)
		blast_record = NCBIXML.read(result_handle)
		counter = 0
		#print("Found ",len(blast_record.alignments), " blast alignments")
		for alignment in blast_record.alignments:
			#print("Retaining alignment for filtering:", alignment.title)
			total_score = 0
			if counter < 1:
				counter = counter + 1
				for hsp in alignment.hsps:
					#print("e value:", hsp.expect, "score:", hsp.score)
					total_score = total_score + hsp.score
				if total_score > 100:
					#print("Retained alignment", alignment.title, "with score", total_score)
					filtered_list.append(seq_record)
SeqIO.write(filtered_list, args['outfile'], "fasta")
