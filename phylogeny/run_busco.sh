#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=2-24:00:00
#SBATCH --mem-per-cpu=12000
#SBATCH --job-name=BUSCO
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 2 ]
then
    	echo -e "\nusage: `basename $0` <reference genome in FASTA format> <BUSCO lineage database>\n"
        echo -e "DESCRIPTION: This script takes a reference genome assembly in FASTA format and estimates the level of genome completeness using BUSCO\n\n"

        echo -e "INPUT:         <reference genome in FASTA format>      Reference genome assembly in FASTA format"
        echo -e "               <BUSCO lineage database>                Path to busco lineage\n\n"

        echo -e "OUTPUT:        <completeness score in BUSCO>           A folder with BUSCO completeness score results\n\n"

        echo -e "REQUIRES:      Activation of BUSCO conda environment\n\n"

        exit
fi


genome=$1
lineage=$2


outDir=$(basename $genome | sed 's/.fasta/.busco/g' |  sed 's/.fa/.busco/g')

echo $outDir
# Run BUSCO with the metaeuk mode
busco -i $genome -o $outDir -m genome -l $lineage --offline -c 10

