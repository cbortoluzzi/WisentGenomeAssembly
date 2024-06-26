#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=FastQC
#SBATCH --output=output_%J
#SBATCH --error=error_%J



if [ $# -ne 3 ]
then
    	echo -e "\nusage: `basename $0` <sample name> <paired-end R1 in FASTQ format> <paired-end R2 in FASTQ format>\n"
        echo -e "DESCRIPTION: This script takes a pair of fastq files and outputs a QC report\n\n"

        echo -e "INPUT:         <sample name>                         Code or name to identify the sample"
        echo -e "               <paired-end R1 in FASTQ format>       First paired-end run in FASTQ format"
        echo -e "               <paired-end R2 in FASTQ format>       Second paired-end run in FASTQ format\n\n"

        echo -e "OUTPUT:        <two reports, one per paired-end run, and one for both paired-end runs combined>\n\n"

	echo -e "REQUIRES:      Requires FastQC (v0.11.9) module available from PATH\n\n"

        exit
fi



sample=$1
fastq_1=$2
fastq_2=$3



# Load module
module load fastqc/0.11.9


mkdir -p fastqc

# This command line runs FastQC, a high throughput sequence QC analysis tool
# The idea here is to produce for each paired-end run in FASTQ format a report that will be used to check for potential problems in the data
fastqc $fastq_1 $fastq_2 -o fastqc --extract -f fastq -d /tmp/

# We will also run fastq-stat as this program outputs some additional statistics that might be of interest
fastq-stat --sampleid $sample --result fastqc/$sample.csv --plot $fastq_1 $fastq_2

