#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=5-24:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=Alignment
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 5 ]
then
    	echo -e "\nusage: `basename $0` <reference genome in FASTA format> <paired-end run R1 in FASTQ format> <paired-end run R2 in FASTQ format> <sample name> <output directory>\n"
        echo -e "DESCRIPTION: This script aligns the paired-end runs of an individual to the reference genome to generate an alignment in BAM format\n\n"

        echo -e "INPUT:         <reference genome in FASTA format>      Reference genome assembly in FASTA format"
        echo -e "               <paired-end run R1 in FASTQ format>     First paired-end run in FASTQ format""
        echo -e "               <paired-end run R2 in FASTQ format>     Second paired-end run in FASTQ format"
        echo -e "               <sample name>                           Name of sample. This is the name that will appear in the VCF file"
	      echo -e "		            <output  directory>		                  Name of the output directory\n\n"

        echo -e "OUTPUT:        <an aligned and sorted BAM file for each paired-end runs>"
        echo -e "               <a QC report>\n\n"

        echo -e "REQUIRES:      Requires BWA (v0.7.17), sambamba (v0.8.1), samblaster (v0.1.24), bamtools (v.2.5.1), java (v14.0.2), and QualiMap (v2.3) available from PATH\n\n"

        exit
fi


export PATH=/cluster/work/pausch/cbortoluzzi/softwares/qualimap_v2.3:$PATH


# Load modules
module load bwa/0.7.17
module load sambamba/0.8.1
module load samblaster/0.1.24
module load bamtools
module load openjdk/14.0.2


ref=$1
fastq_1=$2
fastq_2=$3
sample=$4
outDir=$5


mkdir -p $outDir/bam


name=$(basename $fastq_1 | sed 's/_1.fastq.gz//g' | sed 's/_R1.fastq.gz//g')

flowcell=`zcat $fastq_1 | head -n 1 | cut -f3 -d':'`
lane=`zcat $fastq_1 | head -n 1 | cut -f4 -d':'`



# Make sure that reference genome is unmasked
unmasked=`echo $ref | sed 's/.fasta//g' | sed 's/.fa//g'`
if [[ ! -f $unmasked.rm.fasta ]];then
	echo -e "Unmasked version of the reference genome not found. Making it now\n\n"
	awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' $ref > $unmasked.rm.fasta
fi


# Create index files for reference genome
if [[ ! -f $unmasked.rm.fasta.bwt ]];then
	echo -e "Index database sequences in the FASTA format not found. Making them now\n\n"
	bwa index $unmasked.rm.fasta
fi


# Align paired-end reads to reference genome (the reference genome is unmasked!)
# We will couple BWA mem with samblaster to mark duplicate reads and samtools view to output a new BAM file
if [[ ! -f $outDir/bam/$name.bam ]];then
	echo -e "Aligned BAM file not found for $sample. Making it now...this might take some time\n\n"
	bwa mem -t 16 -T 20 -R "@RG\tID:$flowcell.$lane\tPL:illumina\tSM:$name\tPU:unknown\tCN:ETH" $unmasked.rm.fasta $fastq_1 $fastq_2 | samblaster | samtools view -Sb - > $outDir/bam/$name.bam
fi

# I have lowered the minimum score to output to 20 [default: 30]

# The options in samblaster are the follwing:
#     -e: exclude reads marked as duplicates from discordant, splitter, and/or unmapped file
#     -d: Output discordant read pairs to this file
#     -s: Output split reads to this file abiding by parameters below
#     -u: Output unmapped/clipped reads as FASTQ to this file abiding by parameters below

## Fast coordinate sorting with sambamba
if [[ ! -f $outDir/bam/$name.sort.bam ]];then
	echo -e "Sorted BAM file not found. Making it now\n\n"
	sambamba sort -t 16 --tmpdir /tmp/ -m 10G -o $outDir/bam/$name.sort.bam $outDir/bam/$name.bam
fi


# Get the alignment statistics from the BAM file
if [[ ! -f $outDir/bam/$name.sort.stats.csv ]];then
	echo -e "Statistics not calculated. Calculating them now\n\n"
	bamtools stats -in $outDir/bam/$name.sort.bam > $outDir/bam/$name.sort.stats.csv
fi

# Get quality statistics with Qualimap
mkdir -p $outDir/bam/$name
qualimap bamqc --java-mem-size=16G -bam $outDir/bam/$name.sort.bam -nt 16 -nw 500 -outdir $outDir/bam/$name

echo -e "Done!\n\n"
