#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=Heterozygosity
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 4 ]
then
    	echo -e "\nusage: `basename $0` <VCF> <sample> <bam> <genome coverage>\n"
        echo -e "DESCRIPTION: This script calculates the genome-wide heterozygosity in each individual using a window size of 1 MbF\n\n"

        echo -e "INPUT:         <VCF>                   A VCF input file (all chromosomes)"
        echo -e "               <sample>                Name of the sample to calculate heterozygosity on - as it appears in the VCF file"
        echo -e "               <bam>                   BAM file of the sample analysed"
        echo -e "               <genome coverage>       Genome coverage as estimated by Qualimap or samtools depth\n\n"

        echo -e "OUTPUT:        <A VCF file for the sample analysed>"
        echo -e "               <A text file with heterozygosity calculated in 1 Mb window>\n\n"

        echo -e "REQUIRES:      Requires VCFTools (v0.1.16) available from PATH\n\n"

        exit
fi



# Load modules
module load vcftools
module load bcftools


vcf=$1
sample=$2
bam=$3
coverage=$4


# Obtain a VCF file for the sample of interest
mkdir -p heterozygosity
echo -e "Calculate heterozygosity using a window size of 1 Mb\n\n"
#vcftools --gzvcf $vcf --indv $sample --recode --recode-INFO-all --out heterozygosity/$sample.vcf
#bgzip -c heterozygosity/$sample.vcf.recode.vcf > heterozygosity/$sample.vcf.gz
#tabix -p vcf heterozygosity/$sample.vcf.gz
#rm heterozygosity/*.recode.vcf

# Run python script to calculate heterozygosity using a 1Mb sliding window approach
python3 calculate_genome_heterozygosity.py --vcf heterozygosity/$sample.vcf.gz --bam $bam --d $coverage --w 10000 --o heterozygosity

echo -e "Done\n\n"

