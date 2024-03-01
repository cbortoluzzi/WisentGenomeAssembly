#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 12
#SBATCH --time=5-24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=Pggb
#SBATCH --output=output_%J
#SBATCH --error=error_%J



# Pangenome Sequence Naming : to change the sequence names according to PanSN-spec, we use fastix:
ls input/*.fasta | while read f; do
        sample_name=$(echo $f | cut -f 1 -d '.')
        fastix -p "${sample_name}#1#" $f >> Bovidae.fasta
done
bgzip -@ 4 Bovidae.fasta
samtools faidx Bovidae.fasta.gz


# Mash-based partitioning
mash dist Bovidae.fasta.gz Bovidae.fasta.gz -s 10000 -i > Bovidae.distances.tsv
python3 mash2net.py -m Bovidae.distances.tsv
python3 net2communities.py -e Bovidae.distances.tsv.edges.list.txt -w Bovidae.distances.tsv.edges.weights.txt -n Bovidae.distances.tsv.vertices.id2name.txt

seq 0 27 | while read i; do
        chromosomes=$(cat Bovidae.distances.tsv.edges.weights.txt.community.$i.txt | cut -f 3 -d '#' | sort | uniq | tr '\n' ' ');
        echo "community $i --> $chromosomes";
done


# Data partitioning : each community can be managed by pggb independently of the others. To partition the communities, execute:
seq 0 27 | while read i; do
        echo "community $i"
        samtools faidx Bovidae.fasta.gz $(cat Bovidae.distances.tsv.edges.weights.txt.community.$i.txt) | bgzip -@ 4 -c > Bovidae.community.$i.fa.gz
        samtools faidx Bovidae.community.$i.fa.gz
done


# Run pggb on each community
seq 0 27 | while read i;do
  singularity run pggb_latest.sif pggb -i Bovidae.community.$i.fa.gz -o Bovidae.community.$i.out -t 16 -s 75000 -p 90 -n 8 -k 31 --skip-viz
done

