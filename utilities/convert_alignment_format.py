# Convert alignment from fasta to stockholm format (for rapidnj)
# February 20 2020
# Jacob E. Lemieux

from Bio import AlignIO
import argparse

parser = argparse.ArgumentParser(description = 'Convert alignment from one format (fasta) to another (stockholm)')
parser.add_argument('-i', '--infile', help = 'input file', required=True)
parser.add_argument('-o', '--outfile', help = 'output filename', required=True)
args = vars(parser.parse_args())

AlignIO.convert(args['infile'], "fasta", args['outfile'], "stockholm")
