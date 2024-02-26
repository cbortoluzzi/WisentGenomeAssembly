#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=4-24:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=Admixture
#SBATCH --output=output_%J
#SBATCH --error=error_%J



if [ $# -ne 1 ]
then
                echo -e "\nusage: `basename $0` <VCF>\n"
                echo -e "DESCRIPTION: This script runs the ADMIXTURE analysis to estimate individual ancestries by efficiently computing maximum likelihood estimates in a parametric model\n\n"

                echo -e "INPUT:           <VCF>                 A VCF input file (all chromosomes)\n\n"

                echo -e "OUTPUT:          <Individual ancestries as estimated by ADMIXTURE>\n\n"

                echo -e "REQUIRES:       Requires Plink (v1.90b6.18), Admixture (v1.3.0), and VCFtools (0.1.16) available from PATH\n\n"

                exit
fi


export PATH=/path/to/admixture_linux-1.3.0:$PATH


# Load modules
module load plink
module load vcftools

# The VCF file contains only bi-allelic SNPs
vcf=$1


# Prepare input for ADMIXTURE: the program requires unlinked (i.e. LD pruned) SNPs in plink format
echo -e "Prepare input for ADMIXTURE: remove water buffalo sample and filter bi-allelic SNPs for linkage disequilibrium"
mkdir -p admixture
zcat $vcf | head -n 1000 | grep '^#' | tail -n1 | cut -f 10- | tr '\t' '\n' > admixture/samples.list
cat admixture/samples.list | awk '{if($1 == "Water_buffalo")next}{print $1"\t"$1}' > admixture/samples.keep

# We are going to filter bi-allelic SNPs for LD, a genotype call of 70% (no more than 30% missing data), and a minor allele frequency of 5%
vcftools --gzvcf $vcf --plink --out admixture/input
plink --file admixture/input --make-bed --chr-set 29 --allow-extra-chr --double-id --keep admixture/samples.keep --recode --indep-pairwise 50 10 0.2 --geno 0.3 --maf 0.05 --out admixture/input
plink --bfile admixture/input --double-id --chr-set 29 --allow-extra-chr --threads 10 --exclude admixture/input.prune.out --make-bed --out admixture/input.admixture

for i in $(seq 2 10);do admixture --cv=10 -B1000 admixture/input.admixture.bed $i | tee admixture/log.$i.out;done

# ADMIXTURE produces 2 outputs: .Q which contains cluster assignments for each individual, and .P which containts for each SNP the population allele frequencies
awk '/CV/ {print $3,$4}' admixture/*out | cut -c 4,7-20 > admixture/cv.error
