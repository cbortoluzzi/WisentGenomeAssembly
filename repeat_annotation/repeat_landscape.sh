#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=RepeatMasker
#SBATCH --output=output_%J
#SBATCH --error=error_%J


align=$1
repeat_library=$2


#This creates an additional file "new_alignment_file.align" which contains the added Kimura divergence field after each alignment.
../../softwares/RepeatMasker/util/calcDivergenceFromAlign.pl -s example.divsum -a new_alignment_file.align $align

# Create a Repeat Landscape graph using the divergence summary data generated with the calcDivergenceFromAlign.pl script.
../../softwares/RepeatMasker/util/createRepeatLandscape.pl -div example.divsum -g 2628394923 > example.html
