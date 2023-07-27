#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=5-24:00:00
#SBATCH --mem-per-cpu=40000
#SBATCH --job-name=RepeatModeler
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 1 ]
then
    	echo -e "\nusage: `basename $0` <reference genome in FASTA format>\n"
        echo -e "DESCRIPTION: This script takes a reference genome assembly in FASTA format and builds a de novo repeat library from scratch with RepeatModeler\n\n"

        echo -e "INPUT:         <reference genome in FASTA format>      Reference genome assembly in FASTA format\n\n"

        echo -e "OUTPUT:        <a de novo repeat library>      	A species-specific de novo repeat library\n\n"

        echo -e "REQUIRES:      Requires RepeatModeler (v2.0.4), samtools (v1.15.1) and Genometools (v1.6.2) accessible from PATH\n\n"

        exit
fi




reference=$1
output_directory=$2


name=$(basename $reference | sed 's/.fasta//g' | sed 's/.fa//g' | sed 's/.fasta.gz//g' | sed 's/.fa.gz//g')


mkdir -p repeat_annotation


# Prepare reference-genome.fasta
unmasked=$name.rm.unmasked.fasta
python3 /cluster/work/pausch/cbortoluzzi/bin/repeat_annotation_pipeline/prepare_fasta.py --i $reference --o repeat_annotation


# Build RepeatModeler BLAST database
name=$(basename $unmasked .fasta)
cd repeat_annotation
if [[ ! -f $name.nhr ]]; then
        echo -e "BLAST database not found. Building it...this might take some time\n\n"
        BuildDatabase -name $name -engine ncbi -dir .
fi


# Run RepeatModeler
if [[ ! -f "RM_*/consensi.fa" ]]; then
        echo -e "De novo repeat identification not found. Making it now...this might take a few days\n\n"
        RepeatModeler -threads 16 -database $name
fi

