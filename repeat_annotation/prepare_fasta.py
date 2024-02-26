#!/usr/bin/env python


# Author : cbortoluzzi@ethz.ch


import sys
import os
import gzip
import argparse
from Bio import SeqIO
from pathlib import Path


parser = argparse.ArgumentParser(description = 'Rename header of reference genome assembly and generate unmasked FASTA file')
parser.add_argument('--i', help = 'Reference genome assembly in FASTA format')
parser.add_argument('--o', help = 'Output directory')



def parse_assembly_report(assembly_report):
	mydict = {}
	with open(assembly_report) as f:
		for line in f:
			# Skip first couple of commented lines
			if not line.startswith("#"):
				line = line.strip().split('\t')
				# Save the Genbank accn as key and the chromosome number as value
				if line[1] == "assembled-molecule":
					mydict[line[4]] = line[2]
				else:
					mydict[line[4]] = line[0]
	return mydict



def change_header_naming_fasta(mydict, fasta, output):
	open_f = open_fasta(fasta)
	with open(output, "w") as output_f:
		for record in SeqIO.parse(open_f, "fasta"):
			header = record.id
			try:
				new_header = mydict[header]
				record.id = new_header
				record.description = ''
				#record.seq = record.seq.upper()
				SeqIO.write(record, output_f, "fasta")
			except KeyError:
				record.id = header
				record.description = ''
				#record.seq = record.seq.upper()
				SeqIO.write(record, output_f, "fasta")



def unmask_genome(fasta, output):
	open_f = open_fasta(fasta)
	with open(output, 'w') as output_f:
		for record in SeqIO.parse(open_f, 'fasta'):
			record.description = ''
			#record.seq = record.seq.upper()
			SeqIO.write(record, output_f, "fasta")



def open_fasta(file):
	if file.endswith('.fasta') or file.endswith('.fa'):
		handle = open(file)
	else:
		handle = gzip.open(file, 'rt')
	return handle
	



if __name__ == "__main__":
	args = parser.parse_args()
	p = Path(args.o)
	p.mkdir(parents = True, exist_ok = True)
	output = Path(args.i).stem.replace('.fasta','').replace('.fa', '') + '.rm.unmasked.fasta'
	output_f = Path(args.o, output)
	if output_f.is_file():
		sys.exit("Terminating the code - the output file already exists!")	
	else:
		# Obtain the `assembly_report.txt` file -> this is necessary for changing the header of the FASTA file
		for f in os.listdir(Path(args.i).parents[0]):
			if f.endswith('_assembly_report.txt'):
				assembly_report = Path(Path(args.i).parents[0], f)
		try:
			# If the `assembly_report` file is found, then proceed
			if assembly_report.is_file():
				chromosome_naming = parse_assembly_report(assembly_report)
				change_header = change_header_naming_fasta(chromosome_naming, args.i, output_f)
		# In case the `assembly_report` doesn't exist, prepare a different output file
		# I expect genome assemblies generated in-house to not have an `assembly_report` file
		except NameError:
			# We will generate just an unmasked verion of the input file, because the sequence header is correct
			unmask = unmask_genome(args.i, output_f)

