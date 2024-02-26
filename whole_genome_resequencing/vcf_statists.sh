#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=40000
#SBATCH --job-name=VCFStats
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 1 ]
then
    	echo -e "\nusage: `basename $0` <VCF>\n"
        echo -e "DESCRIPTION: This script calculates a series of statistics on the filtered VCF file\n\n"

        echo -e "INPUT:         <VCF>   A filtered VCF file\n\n"

        echo -e "OUTPUT:        <A set of statistics as calculated by VCFtools>\n\n"

        echo -e "REQUIRES:      Requires VCFTools (v0.1.16) available from PATH\n\n"

        exit
fi




vcf=$1


# Load modules
module load vcftools


mkdir -p vcf_stats


echo -e "Calculate SNP density using a window of 100 Kb\n\n"
vcftools --gzvcf $vcf --SNPdensity 100000 --out vcf_stats/$(basename $vcf)


echo -e "Calculate per site SNP quality (as found in the QUAL column of the VCF file)\n\n"
vcftools --gzvcf $vcf --site-quality --out vcf_stats/$(basename $vcf)


echo -e "Calculte missingness on a per-individual basis\n\n"
vcftools --gzvcf $vcf --missing-indv --out vcf_stats/$(basename $vcf)


echo -e "Calculate missingness on a per-site basis\n\n"
vcftools --gzvcf $vcf --missing-site  --out vcf_stats/$(basename $vcf)


echo -e "Calculate mean depth per individual\n\n"
vcftools --gzvcf $vcf --depth --out vcf_stats/$(basename $vcf)

