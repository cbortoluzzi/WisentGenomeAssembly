#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=40000
#SBATCH --job-name=Picard
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 3 ]
then
    	echo -e "\nusage: `basename $0` <alignment in BAM format> <read group> <sample name>\n"
        echo -e "DESCRIPTION: This script renames the read group of an aligned BAM file\n\n"

        echo -e "INPUT:         <alignment in BAM format>             Alignment in BAM format for which the read group must be changed"
        echo -e "               <read group>                          New read group code/name"
        echo -e "               <sample name>                         New sample code/name\n\n"

        echo -e "OUTPUT:        <two reports, one per paired-end run, and one for both paired-end runs combined>\n\n"
        
       	echo -e "REQUIRES:      Requires Picard (v2.25.7) and Samtools (v1.15.1) module available from PATH\n\n"

        exit
fi


bam=$1
RGID=$2
RGSM=$3


# Load modules
module load picard/2.25.7
module load samtools/1.15.1


# Run AddOrReplaceReadGroups on the alignment in BAM format
# Example:
# picard AddOrReplaceReadGroups I=WISENTM_Urim.sort.bam O=WISENTM_Urim.sort.rm.bam RGID=HWNJ3DSX10.1 RGPL=illumina RGLB=1 RGPU=unknown RGSM=WISENTM_Urim
output=`echo $bam | sed 's/.sort.bam/.sort.rm.bam/g'`
picard AddOrReplaceReadGroups I=$bam O=$output RGID=$RGID RGPL=illumina RGLB=1 RGPU=unknown RGSM=$RGSM
samtools index -@ 10 $output


