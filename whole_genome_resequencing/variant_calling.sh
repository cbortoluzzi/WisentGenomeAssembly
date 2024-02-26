#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=1-24:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=VariantCalling
#SBATCH --output=output_%J
#SBATCH --error=error_%J




if [ $# -ne 3 ]
then
	echo -e "\nusage: `basename $0` <reference genome in FASTA format> <list bam files> <region or chromosome>\n"
	echo -e "DESCRIPTION: This script runs Freebayes for variant calling and a combination of BCFtools and VCFtools for variant postfiltering\n\n"

	echo -e "INPUT:         <reference genome in FASTA format>	Reference genome assembly in FASTA format"
	echo -e "		            <list bam files>			A list of BAM files on which to call variants, one per line"
	echo -e "               <region or chrosome> 			A region (in the form of chr:start-end) or a chromosome\n\n"

	echo -e "OUTPUT:        <A filtered VCF file with called genotypes>"
	echo -e "               <A summary statistics file for each filtered VCF file>\n\n"

	echo -e "REQUIRES:      Requires freebayes (v0.9.21), BCftools (v1.15.1) and VCFtools (v0.1.16) available from PATH\n\n"
  exit
fi



# Load modules
module load bcftools/1.15.1
module load vcftools/0.1.16
module load htslib/1.15.1


ref=$1
bam=$2
reg=$3


mkdir -p vcf


# Variant calling with Freebayes
if [[ ! -f vcf/chromosome.$reg.vcf ]];then
	echo -e "Perform variant calling with Freebayes...this might take some time\n\n"
	freebayes --bam-list $bam --fasta-reference $ref --region $reg --vcf vcf/chromosome.$reg.vcf --ploidy 2 --haplotype-length 0 --min-mapping-quality 20 --min-base-quality 20 --min-alternate-fraction 0.20 --min-alternate-count 2
fi


# Run a custom python script to filter out genotype outliers using the coverage information
if [[ ! -f vcf/chromosome.$reg.fgeno.vcf.gz ]]; then
	echo -e "Filter genotypes based on individual average genome-wide coverage\n\n"
	python3 filter_freebayes_by_coverage.py --vcf vcf/chromosome.$reg.vcf
	bgzip -c vcf/chromosome.$reg.fgeno.vcf > vcf/chromosome.$reg.fgeno.vcf.gz
	tabix -p vcf vcf/chromosome.$reg.fgeno.vcf.gz
fi


# Filter variants based on phred-quality score (>30), allele count (>2), position with respect to indels (5 bases)
if [[ ! -f vcf/chromosome.$reg.fgeno.qual.vcf.gz ]];then
	echo -e "Filter genotypes in BCFtools based on PHRED-quality score, allele count, and position relative to InDels\n\n"
	bcftools filter -i 'QUAL>30 & AC > 2' --IndelGap 5 --SnpGap 5 --threads 10 vcf/chromosome.$reg.fgeno.vcf.gz -Oz -o vcf/chromosome.$reg.fgeno.qual.vcf.gz
	tabix -p vcf vcf/chromosome.$reg.fgeno.qual.vcf.gz
fi


# Generate a file reporting the missingness on a per-site basis and discard sites where all samples have a missing genotype
if [[ ! -f vcf/chromosome.$reg.fgeno.qual.rm.vcf.gz ]];then
	echo -e "Filter sites where all individuals have a missing genotype\n\n"
	vcftools --gzvcf vcf/chromosome.$reg.fgeno.qual.vcf.gz --missing-site --out vcf/chromosome.$reg

	cat vcf/chromosome.$reg.lmiss | sed '1d' | awk '{if($6 == 1)print}' > to_remove.$reg.txt
	vcftools --gzvcf vcf/chromosome.$reg.fgeno.qual.vcf.gz --exclude-positions to_remove.$reg.txt --recode --recode-INFO-all --stdout | bgzip -c > vcf/chromosome.$reg.fgeno.qual.rm.vcf.gz
	tabix -p vcf vcf/chromosome.$reg.fgeno.qual.rm.vcf.gz
fi


# Calculate statistics
if [[ ! -f vcf/chromosome.$reg.fgeno.qual.rm.stats ]];then
	echo -e "Calculate statistics on filtered VCF file\n\n"
	bcftools stats vcf/chromosome.$reg.fgeno.qual.rm.vcf.gz --threads 10 > vcf/chromosome.$reg.fgeno.qual.rm.stats
fi


echo -e "Done\n\n"



