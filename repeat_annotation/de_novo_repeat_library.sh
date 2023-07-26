#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=5-24:00:00
#SBATCH --mem-per-cpu=40000
#SBATCH --job-name=RepeatModeler
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 2 ]
then
	echo -e "\nusage: `basename $0` <reference genome in FASTA format> <output directory>\n"
	echo -e "DESCRIPTION: This script takes a reference genome assembly in FASTA format and builds a de novo repeat library from scratch with RepeatModeler\n\n"

	echo -e "INPUT:          <reference genome in FASTA format>      Reference genome assembly in FASTA format"
	echo -e "                <output directory>                      Name of output directory\n\n"
 	
  	echo -e "OUTPUT:        <a RepeatModeler BLAST database>"
	echo -e "               <a de novo repeat library>"

	echo -e "REQUIRES:      Requires RepeatModeler (v2.0.4), samtools (v1.15.1) and Genometools (v1.6.2) accessible from PATH\n\n"

	exit
fi


reference=$1
output_directory=$2


species_name=$(dirname $reference | rev | cut -f 1 -d '/' | rev)
outDir=`echo $output_directory'/'$species_name`


mkdir -p $outDir


# Prepare reference-genome.fasta
unmasked=$(basename $reference | sed 's/.fasta.gz//g' | sed 's/.fa.gz//g')".rm.unmasked.fasta"
python3 /cluster/work/pausch/cbortoluzzi/bin/repeat_annotation_pipeline/prepare_fasta.py --i $reference --o $outDir


# Build RepeatModeler BLAST database
name=$(basename $unmasked .fasta)
mkdir -p $outDir/RepeatModeler
cd $outDir/RepeatModeler
if [[ ! -f $name.nhr ]]; then
	echo -e "BLAST database not found. Building it...this might take some time\n\n"
	BuildDatabase -name $name -engine ncbi -dir $outDir
fi


# Run RepeatModeler
if [[ ! -f "RM_*/consensi.fa" ]]; then
	echo -e "De novo repeat identification not found. Making it now...this might take a few days\n\n"
	RepeatModeler -threads 16 -database $name
fi

