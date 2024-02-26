#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=1-24:00:00
#SBATCH --mem-per-cpu=12000
#SBATCH --job-name=compleasm
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 2 ]
then
    	echo -e "\nusage: `basename $0` <reference genome in FASTA format> <BUSCO lineage database>\n"
        echo -e "DESCRIPTION: This script takes a reference genome assembly in FASTA format and estimates the level of genome completeness using compleasm, which relies on miniprot\n\n"

        echo -e "INPUT:         <reference genome in FASTA format>        Reference genome assembly in FASTA format"
        echo -e "               <species name>                            Name of species to run compleasm on\n\n"

        echo -e "OUTPUT:        <completeness score in compleasm>         A folder with compleasm completeness score results\n\n"

        echo -e "REQUIRES:      Activation of compleasm conda environment\n\n"

        exit
fi


genome=$1
species=$2


mkdir -p compleasm/$species
mkdir -p library_path

# Download library
compleasm download -L library_path mammalia

# Run miniprot
miniprot --trans -u -I --outs=0.95 --gff -t 10 $genome library_path/mammalia_odb10/refseq_db.faa.gz > compleasm/$species/output.gff

# Analyze miniprot output
compleasm analyze -g compleasm/$species/output.gff -l mammalia -o compleasm/$species -t 10 -L library_path/ -m busco

