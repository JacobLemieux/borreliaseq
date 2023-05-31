from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re


# to do: modify regular expression to match with/without case and with/without dash

parser = argparse.ArgumentParser(description = 'Extract sequences of interest')
parser.add_argument('-i', '--infile', help = 'input file', required=True)
parser.add_argument('-o', '--outfile', help = 'output filename', required=True)
parser.add_argument('-t', '--filetype', help = 'file type (genbank or fasta)', required =True)
args = vars(parser.parse_args())

nrecords = len(list(SeqIO.parse(args['infile'], args['filetype'])))
seqlist = []

pattern = re.compile("PF-32")
pattern2 = re.compile("PF32")
pattern3 = re.compile("ParA family protein")

plasmid = ""
strain = ""

open(args['outfile'], "w")
if args['filetype'] == "genbank":
	for seq_record in SeqIO.parse(args['infile'], "genbank"):
		#print(len(seq_record.features))
		for feature in seq_record.features:
			if feature.type == "source" and ("plasmid" in feature.qualifiers.keys()):
				#print(feature)
				plasmid = feature.qualifiers["plasmid"][0]
				plasmid = plasmid.replace(" ", "_")
				if "strain" in feature.qualifiers.keys():
					strain = feature.qualifiers["strain"][0]
			if feature.type == "CDS" and ("translation" in feature.qualifiers.keys()):  # this approach is error-prone because it assumes source annotated as a plasmid is always followed by partition protein
				#print(str(feature.qualifiers["product"]))
				if pattern.search(str(feature.qualifiers["product"])) or pattern2.search(str(feature.qualifiers["product"])) or pattern3.search(str(feature.qualifiers["product"])):
					#print(feature.extract(seq_record))
					loc_name = feature.qualifiers["locus_tag"][0].split("_")[1]
					new_seq = SeqRecord(Seq(feature.qualifiers["translation"][0]), id = loc_name + "_"+ strain + "_" + plasmid, description = "") # + "_"+str(len(seq_record.seq)) 
					#print(feature.qualifiers)
					if len(seq_record.seq) < 150000:
						seqlist.append(new_seq)
if args['filetype'] == "fasta":
	for seq_record in SeqIO.parse(args['infile'], "fasta"):	
		if pattern.search(seq_record.description) or pattern2.search(seq_record.description) or pattern3.search(seq_record.description):
			#print(seq_record.description)
			loc_name = seq_record.description.split(" ")[0].split("_")[1]
			new_seq = SeqRecord(seq_record.seq, id = args['infile'].split("/")[1] + "_" + loc_name, description = "")
			seqlist.append(new_seq)
SeqIO.write(seqlist, args['outfile'], "fasta")