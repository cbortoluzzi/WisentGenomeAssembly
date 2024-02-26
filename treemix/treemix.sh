#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=3-24:00:00
#SBATCH --mem-per-cpu=12000
#SBATCH --job-name=Treemix
#SBATCH --output=output_%J
#SBATCH --error=error_%J


vcf=$1


module load plink
module load bcftools
module load vcftools


# We are going to run Treemix. Treemix assumes unlinked SNPs, so we are first going to prune the input VCF file for SNPs in high LD (--indep-pairwise option).
# Treemix does not like missing data, so we will remove sites with missing data as well (--geno option)
echo -e "Prune linked SNPs\n\n"
vcftools --gzvcf $vcf --plink --out input
plink --file input --make-bed --chr-set 29 --allow-extra-chr --double-id --indep-pairwise 50 10 0.2 --geno 0 --out input
plink --bfile input --make-bed --chr-set 29 --allow-extra-chr --double-id --exclude input.prune.out --out input.treemix

# Generate cluster file
echo -e "Generate cluster sample file\n\n"
bcftools query -l $vcf |  awk '{split($1,pop,"."); print $1"\t"$1"\t"$1}' > samples.clust

# Create input for Treemix
echo -e "Create input file for Treemix\n\n"
plink --bfile input.treemix --freq --missing --within samples.clust --out input.treemix --chr-set 29 -allow-extra-chr
gzip input.treemix.frq.strat

python2.7 plink2treemix.py input.treemix.frq.strat.gz input.treemix.frq.gz

# Run Treemix
echo -e "Run Treemix\n\n"
for i in $(seq 1 5);do treemix -i input.treemix.frq.gz -o input.treemix.frq.$i -k 1000 -m $i -root Water_buffalo -se -bootstrap 1000 -noss > treemix\_${i}\_log; done

