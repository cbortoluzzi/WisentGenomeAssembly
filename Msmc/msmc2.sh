#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=2-24:00:00
#SBATCH --mem-per-cpu=12000
#SBATCH --job-name=MSMC2
#SBATCH --output=output_%J
#SBATCH --error=error_%J



if [ $# -ne 3 ]
then
    	echo -e "\nusage: `basename $0` <reference genome in FASTA format> <bam> <genome coverage>\n"
      echo -e "DESCRIPTION: This script estimates the demographic history of an individual using the MSMC2 approach\n\n"

      echo -e "INPUT:         <reference genome in FASTA format>  A reference genome in FASTA format"
      echo -e "               <bam>                               A bam file of the sample analysed"
      echo -e "               <genome coverage>                   Genome coverage as estimated by Qualimap or samtools depth\n\n"

      echo -e "OUTPUT:        <A reconstructed demographic history>\n\n"

      echo -e "REQUIRES:      Requires BCFtools (1.16), bwa (v0.7.17), msmc2, and R (v4.2.2) available from PATH\n\n"

      exit
fi


export PATH=/path/to/msmc2/build/release:$PATH
export PATH=/path/to/seqbility-20091110:$PATH


ref=$1
bam=$2
coverage=$3


# Load module
module load r/4.2.2
module load bcftools
module load bwa/0.7.17


# Generate consensus sequences and mask files for the sample of interest
echo -e "Generate consensus sequences and mask file for the sample analysed\n\n"
sample=$(basename $bam | sed 's/.sort.bam//g' | sed 's/.sort.rm.bam//g')
for i in $(seq 1 29);do
	# By default in bamCaller, the minMapQ is 20 and the minConsQ is 20
	bcftools mpileup -q 20 -Q 20 -C 50 --threads 10 -r $i -f $ref $bam | bcftools call -c -V indels | /path/to/msmc2/msmc-tools/bamCaller.py $coverage $sample.$i.mask.bed.gz | gzip -c > $sample.$i.vcf.gz
	/path/to/msmc2/msmc-tools/generate_multihetsep.py --chr $i --mask $sample.$i.mask.bed.gz $sample.$i.vcf.gz > $sample.$i.multihetsep.txt
done

# Combine all inputs in one single file
list=`for i in $(seq 1 29);do echo $sample.$i.multihetsep.txt;done`
# We will use 50 expectation maximization iteration (-i 50)
echo -e "Run MSMC2\n\n"
msmc2 -i 50 -t 10 -p 1*2+25*1+1*2+1*3 -o $sample $list


# Plot demographic history
Rscript plot_msmc2.r $sample.final.txt


