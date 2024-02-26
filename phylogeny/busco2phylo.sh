#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=5-24:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=busco2phylo
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 1 ]
then
	echo -e "\nusage: `basename $0` <list of assemblies>\n"
	echo -e "DESCRIPTION: This script builds a phylogenetic tree from a set of complete, single copy BUSCO genes\n\n"

	echo -e "INPUT:         <list of assemblies>     A text file with the name of all assemblies to include in the phylogenetic tree, one per line\n\n"
 
	echo -e "OUTPUT:        <a supermatrix>          A supermatrix containing the aa sequence of all complete, single copy BUSCO genes\n"
 	echo -e "               <a phylogenetic tree>    A phylogenetic tree obtained from RAxML after 1,000 bootstraps\n\n"

	echo -e "REQUIRES:      Requires mafft (v7.490), trimAl (v1.4.rev22), and RAxML (v8.2.12) available from PATH\n\n"

  exit
fi


# Export paths
export PATH=/path/to/mafft-7.490-with-extensions/bin:$PATH
export PATH=/path/to/trimal/source:$PATH
export PATH=/path/to/standard-RAxML:$PATH



list_assembly=$1


num_genomes=`wc -l $list_assembly | awk '{print $1}'`

cat $list_assembly | sed 's/.fasta/.busco/g' |  sed 's/.fa/.busco/g' | while read assembly;do cat $assembly/mammalia_odb10/run_mammalia_odb10/full_table.tsv | grep -v '^#' | awk '$2 == "Complete"{print $1}' >> complete_busco_ids.txt; done
sort complete_busco_ids.txt | uniq -c | awk '$1=="'$num_genomes'"{print $2}' > final_busco_ids.txt


mkdir -p mafft && mkdir -p trimal

cat final_busco_ids.txt | while read busco_id;do

 	# Obtain protein sequence for those BUSCO genes that are comple, sigle copy and are present in all species included in the dataset
	mkdir -p fasta/$busco_id
	cat $list_assembly | sed 's/.fasta/.busco/g' |  sed 's/.fa/.busco/g' | while read assembly;do
		for faa in $assembly/mammalia_odb10/run_mammalia_odb10/busco_sequences/single_copy_busco_sequences/$busco_id.faa
		do
			cat $faa | awk '/^>/{print ">'$species'"; next}{print}' > $outDir/fasta/$busco_id/$species.$busco_id.fa
		done
	cat fasta/$busco_id/*.fa >> mafft/$busco_id.aln

 	# Align protein sequences with MAFFT
	mafft --amino mafft/$busco_id.aln > mafft/$busco_id.aln.mafft
	rm mafft/$busco_id.aln

 	# Trim alignments with trimAL
	trimal -in mafft/$busco_id.aln.mafft -out trimal/$busco_id.aln.mafft.trimal -gt 0.8 -st 0.001
done


# Generate supermatrix
mkdir -p matrix
python3 superalignment.py --i trimal --o matrix


# Run Raxml with 1,000 bootstraps
raxmlHPC-SSE3 -T 16 -f a -m PROTGAMMAJTT -N 1000 -n my_busco_phylo -s matrix/supermatrix.aln.mafft.trimal.fa -p 13432 -x 89090


