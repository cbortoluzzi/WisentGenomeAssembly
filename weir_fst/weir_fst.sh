#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=Fst
#SBATCH --output=output_%J
#SBATCH --error=error_%J

if [ $# -ne 1 ]
then
    	echo -e "\nusage: `basename $0` <VCF>\n"
        echo -e "DESCRIPTION: This script calculates the Fst statistics between two populations using a window size of 50 Kb\n\n"

        echo -e "INPUT:           <VCF>        A VCF input file (all chromosomes)\n\n"

        echo -e "OUTPUT:          <Fst statistics between two populations>\n\n"

        echo -e "REQUIRES:       Requires VCFTools (v0.1.16) available from PATH\n\n"

        exit
fi



# Load modules
module load vcftools


vcf=$1


# Fst is the proportion of the total genetic variance contained in a subpopulation (the S subscript) relative to the total genetic variance (the T subscript). Values can range from 0 to 1. 
# High FST implies a considerable degree of differentiation among populations.

mkdir -p Fst

# The input VCF file has only bi-allelic SNPs
zcat $vcf | head -n 1000 | grep '^#' | tail -n1 | cut -f 10- | tr '\t' '\n' > Fst/samples.list
# Select European bison samples
cat Fst/samples.list | grep -e 'WISENT' -e 'BBO' > Fst/bison_bonasus.list
# Select American bison samples
cat Fst/samples.list | grep -e 'American_bison' > Fst/bison_bison.list
# Select cattle samples
cat Fst/samples.list | awk '{if($1 == "Angus" || $1 == "Hereford" || $1 == "Charolais" || $1 == "Holstein")print}' > Fst/bos_taurus.list

# Calculate Fst using a window size of 50 Kb
# Bison bonasus vs Bison bison
vcftools --gzvcf $vcf --weir-fst-pop Fst/bison_bonasus.list --weir-fst-pop Fst/bison_bison.list --fst-window-size 50000 --out Fst/Fst.bison_bonasus.vs.bison_bison

# Bison bonasus vs Bos taurus
vcftools --gzvcf $vcf --weir-fst-pop Fst/bison_bonasus.list --weir-fst-pop Fst/bos_taurus.list --fst-window-size 50000 --out Fst/Fst.bison_bonasus.vs.bos_taurus

# Bison bison vs Bos taurus
vcftools --gzvcf $vcf --weir-fst-pop Fst/bison_bison.list --weir-fst-pop Fst/bos_taurus.list --fst-window-size 50000 --out Fst/Fst.bison_bison.vs.bos_taurus
