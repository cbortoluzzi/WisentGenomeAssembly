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
    	echo -e "\nusage: `basename $0` <a de novo repeat library> <reference genome in FASTA format>\n"
        echo -e "DESCRIPTION: This script runs RepeatMasker using as input the de-novo repeat library built from LTR retriever and RepeatModeler\n\n"

        echo -e "INPUT:         <a de novo repeat library>            De-novo repeat library built from RepeatModeler"
        echo -e "               <reference genome in FASTA format>    Reference genome assembly in FASTA format\n\n"


        echo -e "OUTPUT:        <a masked reference genome assembly in FASTA format>\n\n"

        echo -e "REQUIRES:      Requires RepeatMasker (v4.1.5) accessible from PATH\n\n"

        exit
fi



export PATH=/path/to/RepeatMasker:$PATH


library=$1
reference=$2


directory=$(dirname $library)


# Run RepeatMasker
RepeatMasker -e rmblast -pa 10 -s -lib $library -dir $directory -a -xsmall -gff $reference

