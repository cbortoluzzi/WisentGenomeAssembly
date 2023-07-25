#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=FastQC
#SBATCH --output=output_%J
#SBATCH --error=error_%J



if [ $# -ne 4 ]
then
    	echo -e "\nusage: `basename $0` <sample_name> <fastq_1> <fastq_2> <output_directory>\n"
        echo -e "DESCRIPTION: This script takes a pair of fastq files and outputs a QC report\n\n"

        echo -e "INPUT:         <sample_name>           Code or name to identify the sample"
        echo -e "               <fastq_1>               First paired end FASTQ file"
        echo -e "               <fastq_2>              	Second paired end FASTQ file"
        echo -e "               <output_directory>      Name of output directory\n\n"

        echo -e "OUTPUT:        <a report, one per pair of fastq files> The report will be stored in the output directory in a directory defined by the sample name\n\n"

        echo -e "REQUIRES:      Loading of FastQC (v0.11.9) module\n\n"

        exit
fi



sample=$1
fastq_1=$2
fastq_2=$3
outDir=$4


# load module
module load fastqc/0.11.9


mkdir -p $outDir/fastqc

# The first thing we need to do is to run FastQC - a high throughput sequence QC analysis tool
# The idea here is to produce for each fastq file a quality control report that will be used to check for potential problems in the data
fastqc $fastq_1 $fastq_2 -o $outDir/fastqc --extract -f fastq -d /tmp/

# We will also run fastq-stat as this program outputs some interesting statistics
fastq-stat --sampleid $sample --result $outDir/fastqc/$sample.csv --plot $fastq_1 $fastq_2



