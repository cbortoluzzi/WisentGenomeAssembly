#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=5-24:00:00
#SBATCH --mem-per-cpu=40000
#SBATCH --job-name=RepeatMasker
#SBATCH --output=output_%J
#SBATCH --error=error_%J




if [ $# -ne 2 ]
then
    	echo -e "\nusage: `basename $0` <de-novo_repeat_library.fasta> <reference-genome.fasta>\n"
        echo -e "DESCRIPTION: This script runs RepeatMasker using as input the de-novo repeat library built from LTR retriever and RepeatModeler\n\n"

        echo -e "INPUT:         <de-novo_repeat_library.fasta>  De-novo repeat library built from LTR retriever and RepeatModeler\n"
        echo -e "INPUT:         <reference-genome.fasta>        Reference genome assembly in FASTA format\n\n"

        echo -e "OUTPUT:        <a masked reference genome assembly in FASTA format>\n\n"

        echo -e "REQUIRES:      Installation of RepeatMasker (v4.1.5) accessible from PATH\n\n"

        exit
fi


# Export path to RepeatMasker
export PATH=/cluster/work/pausch/cbortoluzzi/softwares/RepeatMasker:$PATH


de_novo_repeat_library=$1
reference_genome=$2


# Run RepeatMasker
outDir=$(dirname $de_novo_repeat_library)
mkdir -p $outDir/RepeatMasker
RepeatMasker -e rmblast -pa 10 -s -lib $de_novo_repeat_library -dir $outDir/RepeatMasker -a -xsmall -gff $reference_genome
