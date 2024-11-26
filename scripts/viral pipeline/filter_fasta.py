# Python
# To get high-quality viruses
# By Cassie

import sys
import Bio
from Bio import SeqIO
import argparse 

parser = argparse.ArgumentParser(description='Process to obtain high-quality genomes from combined fasta file. Created to get viral genome collection of only high-quality genomes from a file containing all predicted viral genomes. Requires file with list of genome names.')
parser.add_argument('-i', '--input', help='A fasta file with multiple sequences')
parser.add_argument('-g', '--genomes', help='A file of genome names')
parser.add_argument('-o', '--out', help='A fasta file name to save output too')

	
args = parser.parse_args()

def filter_fasta(input_fasta, output_file, genomes): 
											
	seq_records = SeqIO.parse(input_fasta, format='fasta') #parses the fasta file
	
	with open(genomes) as f:
			genomes_ids_list = f.read().splitlines() #parse the contamination file which is each line as a scaffold id 
	
	OutputFile = open(output_file, 'w') #opens new file to write to
	
	for record in seq_records: 
		#print(record.id)
		if record.id in genomes_ids_list:
			OutputFile.write('>'+ record.id +'\n') #writes the scaffold to the file (or assession) 
			OutputFile.write(str(record.seq)+'\n') #writes the seq to the file
	
	OutputFile.close()


#Run

filter_fasta(args.input, args.out, args.genomes)
