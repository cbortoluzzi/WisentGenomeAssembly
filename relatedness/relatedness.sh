
#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=relatedness
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 1 ]
then
    	echo -e "\nusage: `basename $0` <VCF>\n"
        echo -e "DESCRIPTION: This script calculates a relatedness statistic (e.g., KING) among all individuals included in the VCF file\n\n"

        echo -e "INPUT:           <VCF>                           A VCF input file (all chromosomes)\n\n"

        echo -e "OUTPUT:          <A filtered VCF file>           A filtered VCF file (bi-allelic SNPs, indels excluded)"
        echo -e "                 <A relatedness statistics>      A relatedness statistic for each individual\n\n"

        echo -e "REQUIRES:        Requires VCFtools (v0.1.16) and BCFtools (v1.6) available from PATH\n\n"

        exit
fi


# Load modules
module load bcftools/1.6
module load vcftools


vcf=$1
output=`echo $vcf | sed 's/.vcf.gz/.biallelic.noindels.vcf.gz/g'`


# Retain only biallelic SNPs and remove indels = this VCF file will be used for many downstream analyses and for further data pruning depending on the type of analysis
echo -e "Obtain a VCF file with only bi-allelic SNPs\n\n"
vcftools --gzvcf $vcf --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > $output


# Calculate relatedness with the --relatedness2 option in VCFTools
echo -e "Calculate relatedness statistic\n\n"
mkdir -p relatedness
vcftools --gzvcf $output --relatedness2 --out relatedness/relatedness

# Plot relatedness as matrix
Rscript plot_relatedness.R relatedness/relatedness.relatedness2




