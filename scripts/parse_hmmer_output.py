### parse hmmer output
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

parser = argparse.ArgumentParser(description = 'Extract sequences of interest')
parser.add_argument('-i', '--infile', help = 'input file', required=True)
parser.add_argument('-m', '--inhmm', help = 'hmmer output that lists hmm hits and is provided as input')
parser.add_argument('-o', '--outfile', help = 'output filename', required=True)
args = vars(parser.parse_args())

hmm = pd.read_csv(args["inhmm"], sep = "\t", comment = "#", header =None)
#hmm = pd.read_csv("/Users/jy17/Dropbox/Tick/borrelia/borreliaseq/results/annotation/short_read2/UMA8_pfam.out", sep = "\t", comment = "#", header =None)
seqs_to_extract = hmm[0].tolist()
seqs_to_extract = [i.split()[0] for i in seqs_to_extract]

outlist = []
for seq_record in SeqIO.parse(args["infile"], "fasta"):
  if seq_record.id in seqs_to_extract:
    outlist.append(seq_record)
    #print(seq_record.id)

SeqIO.write(outlist,args["outfile"], format = "fasta")
