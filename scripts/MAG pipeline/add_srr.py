# Python
# To add SRR numbers and renumber contigs
# By Cassie

import sys
import Bio
from Bio import SeqIO
import argparse 

parser = argparse.ArgumentParser(description='Process to rename contigs and add SRR numbers for deposition to NCBI.')
parser.add_argument('-i', '--input', help='A fasta file with multiple sequences')
parser.add_argument('-s', '--SRR', help='SRR numbers for NCBI')
parser.add_argument('-o', '--out', help='A fasta file name to save output too')

	
args = parser.parse_args()

def add_SRR(input_fasta, SRR_values:str, output_file): 
											
    seq_records = SeqIO.parse(input_fasta, format='fasta') #parses the fasta file
	
    OutputFile = open(output_file, 'w') #opens new file to write to
	
    i = 1

    for record in seq_records: 
		#print(record.id)
        OutputFile.write(f'>contig{i} [SRA={SRR_values}]\n') #writes the scaffold to the file (or assession) 
        OutputFile.write(f'{record.seq}\n') #writes the seq to the file
        i += 1
	
    OutputFile.close()


#Run

add_SRR(args.input, args.SRR, args.out)
