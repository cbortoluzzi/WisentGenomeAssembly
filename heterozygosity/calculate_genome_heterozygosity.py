#!/usr/bin/env python



# Author : @cb46



import vcf
import argparse
import subprocess
from pathlib import Path



parser = argparse.ArgumentParser(description = 'Calculate genome-wide heterozygosity using a sliding window approach')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--bam', help = 'BAM file')
parser.add_argument('--d', help = 'Genome-wide depth', type = float)
parser.add_argument('--w', help = 'Window size [default = 10000 bp]', type = int, default = 10000)
parser.add_argument('--o', help = 'Output directory')



class Heterozygosity:

	def sequences_bam(self, bam_f):
		self.mygenome = {}
		# Get BAM index stats: we will retain only the chromosome and total sequence length (in bp)
		command = 'samtools idxstats %s | cut -f 1,2' %(bam_f)
		cmd = subprocess.check_output(command, shell = True).decode()
		outcmd = cmd.split('\n')
		for line in outcmd:
			if line:
				chromosome, length = line.strip().split()
				length = int(length)
				try:
					if isinstance(int(chromosome), int):
						self.mygenome[chromosome] = length
				except ValueError:
					continue
		return self.mygenome


	def calculate_binned_heterozygosity(self, window, min_depth, max_depth, bam_f, vcf_f, filename, path):
		for chromosome in self.mygenome:
			seq_length = self.mygenome[chromosome]
			for i in range(0, seq_length, window):
				start = i
				end = i + window
				self.bam_depth(chromosome, start, end, min_depth, max_depth, bam_f, vcf_f, window, filename, path)


	def bam_depth(self, chromosome, start, end, min_depth, max_depth, bam_f, vcf_f, window, filename, path):
		cov_sites = 0
		# Obtain the read depth of each site in the BAM file
		command = 'samtools depth -r %s:%d-%d %s' %(chromosome, start, end, bam_f)
		cmd = subprocess.check_output(command, shell = True).decode()
		outcmd = cmd.split('\n')
		for line in outcmd:
			if line:
				chromosome, position, depth = line.strip().split()
				# Filter sites based on read depth
				if int(depth) >= min_depth and int(depth) <= max_depth:
					cov_sites += 1
		self.heterozygosity(chromosome, start, end, cov_sites, vcf_f, window, filename, path)


	def heterozygosity(self, chromosome, start, end, cov_sites, vcf_f, window, filename, path):
		vcf_reader = vcf.Reader(filename=vcf_f)
		nhet = 0
		for record in vcf_reader.fetch(chromosome, start, end):
			for call in record.samples:
				nhet += record.num_het
		if cov_sites == window + 1:
			cov_sites = window
		try:
			SNPcount = round((window / cov_sites) * nhet, 3)
		except ZeroDivisionError:
			SNPcount = 0.0
		output_f = Path(path, filename + '.heterozygosity.txt')
		with open(output_f, 'a') as out:
			out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome, start, end, cov_sites, nhet, SNPcount))




if __name__ == "__main__":
	args = parser.parse_args()
	filename = Path(args.vcf).stem.replace('.vcf', '')
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	genome_wide_heterozygosity = Heterozygosity()
	genome_wide_heterozygosity.sequences_bam(args.bam)
	min_depth = round((1/3) * args.d, 2)
	max_depth = round(2.5 * args.d, 2)
	print ("Minimum depth:", min_depth)
	print ("Maximum depth:", max_depth)
	genome_wide_heterozygosity.calculate_binned_heterozygosity(args.w, min_depth, max_depth, args.bam, args.vcf, filename, args.o)

