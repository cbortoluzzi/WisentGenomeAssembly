#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 20
#SBATCH --time=10-24:00:00
#SBATCH --mem-per-cpu=80000
#SBATCH --job-name=Pggb
#SBATCH --output=output_%J
#SBATCH --error=error_%J


module load bcftools

input=$1 # PGGB graph output

vg deconstruct -P input/Bison_bonasus -a -e -d 1 -t 16 $input | bgzip -c -@ 16 $input.raw.vcf > $input.raw.vcf.gz

bcftools annotate --threads 34 -x INFO $input.raw.vcf.gz | vcfwave -t 34 -L 1000000 -k > $input.wave.vcf

cat $input.wave.vcf | sed 's/input\/Bison_bonasus#1#//g' | bcftools norm --threads 16 -m -any -f input/Bison_bonasus.chr.fasta | bcftools norm --threads 16 -d none | bcftools sort -T /tmp/ -o $input.wave.norm.vcf
bcftools view -i 'abs(ILEN)>=50' $input.wave.norm.vcf > $input.wave.norm.SV.vcf


