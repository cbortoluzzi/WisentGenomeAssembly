#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=5-24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=RepeatModeler
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 2 ]
then
        echo -e "\nusage: `basename $0` <reference genome in FASTA format>\n"
        echo -e "DESCRIPTION: This script takes a reference genome assembly in FASTA format and builds a de novo repeat library from scratch with RepeatModeler\n\n"

        echo -e "INPUT:         <reference genome in FASTA format>      Reference genome assembly in FASTA format"
        echo -e "               <species name>                          Name of species whose reference genome is being analysed\n\n"

        echo -e "OUTPUT:        <a de novo repeat library>      A species-specific de novo repeat library\n\n"

        echo -e "REQUIRES:      Requires RepeatModeler (v2.0.4), samtools (v1.19.2) and Genometools (v1.6.2) accessible from PATH\n\n"

        exit
fi




reference=$1
species=$2


name=$(basename $reference | sed 's/.fasta.gz//g' | sed 's/.fa.gz//g')


mkdir -p repeat_annotation/$species


# Prepare reference-genome.fasta
unmasked=$name.rm.unmasked.fasta
echo -e "Prepare reference genome in FASTA format\n\n"
python3 prepare_fasta.py --i $reference --o repeat_annotation/$species


# Build RepeatModeler BLAST database
database=$(basename $unmasked .fasta)
cd repeat_annotation/$species
if [[ ! -f $database.nhr ]]; then
        echo -e "BLAST database not found. Building it ... this might take some time\n\n"
        BuildDatabase -name $database -engine ncbi -dir .
fi


# Run RepeatModeler
if [[ ! -f "RM_*/consensi.fa" ]]; then
        echo -e "De novo repeat identification not found. Making it now...this might take a few days\n\n"
        RepeatModeler -threads 16 -database $database
fi


