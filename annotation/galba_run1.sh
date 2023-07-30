#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 8
#SBATCH --time=5-24:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=GALBA
#SBATCH --output=output_%J
#SBATCH --error=error_%J



export PATH=/cluster/work/pausch/cbortoluzzi/softwares/GALBA-1.0.7:$PATH
export PATH=/cluster/work/pausch/cbortoluzzi/softwares/GALBA-1.0.7/scripts:$PATH
export PATH=/cluster/work/pausch/cbortoluzzi/softwares/BRAKER-3.0.3/scripts:$PATH
export PATH=/cluster/work/pausch/cbortoluzzi/softwares/Augustus-3.5.0/bin:$PATH
export ALIGNMENT_TOOL_PATH=/cluster/work/pausch/cbortoluzzi/softwares/gth-1.7.3-Linux_x86_64-64bit/bin:$PATH
export AUGUSTUS_CONFIG_PATH=/cluster/work/pausch/cbortoluzzi/softwares/Augustus-3.5.0/config/
export AUGUSTUS_BIN_PATH=/cluster/work/pausch/cbortoluzzi/softwares/Augustus-3.5.0/bin
export AUGUSTUS_SCRIPTS_PATH=/cluster/work/pausch/cbortoluzzi/softwares/Augustus-3.5.0/scripts
export GENEMARK_PATH=/cluster/work/pausch/cbortoluzzi/softwares/GeneMark-ETP/bin/gmes
export TSEBRA_PATH=/cluster/work/pausch/cbortoluzzi/softwares/TSEBRA-v.1.1.1/bin
export CDBTOOLS_PATH=/cluster/work/pausch/cbortoluzzi/softwares/cdbfasta
export PROTHINT_PATH=/cluster/work/pausch/cbortoluzzi/softwares/ProtHint/bin
export BLAST_PATH=/cluster/work/pausch/cbortoluzzi/softwares/ncbi-blast-2.14.0+/bin
export MINIPROT_PATH=/cluster/work/pausch/cbortoluzzi/softwares/miniprot-0.11
export SCORER_PATH=/cluster/work/pausch/cbortoluzzi/softwares/miniprot-boundary-scorer
export MINIPROTHINT_PATH=/cluster/work/pausch/cbortoluzzi/softwares/miniprothint


# Load modules
module load bcftools
module load htslib
module load samtools
module load diamond
module load exonerate
module load bamtools


# Run GALBA using the protein sequence of Bos taurus
echo -e "Run GALBA ... this might take a few days\n\n"
perl /cluster/work/pausch/cbortoluzzi/softwares/GALBA-1.0.7/scripts/galba.pl --species=Bison_bonasus --genome=/cluster/work/pausch/cbortoluzzi/wisent_project/repeat_annotation/Bison_bonasus/RepeatMasker/WISENTM_Urano_F1.hap1.rm.unmasked.fasta.masked --prot_seq=/cluster/work/pausch/alex/REF_DATA/Bos_taurus.ARS-UCD1.2.109.proteins.fa --threads=4 --workingdir=/cluster/work/pausch/cbortoluzzi/wisent_project/annotation/WISENTM_Urano_F1_hap1/run_1

# Reducing noise in augustus.hints.gtf with TSEBRA
echo -e "Reduce noise with TSEBRA\n\n"
python3 /cluster/work/pausch/cbortoluzzi/softwares/TSEBRA-v.1.1.1/bin/tsebra.py -g /cluster/work/pausch/cbortoluzzi/wisent_project/annotation/WISENTM_Urano_F1_hap1/run_1/augustus.hints.gff -e /cluster/work/pausch/cbortoluzzi/wisent_project/annotation/WISENTM_Urano_F1_hap1/run_1/hintsfile.gff -o /cluster/work/pausch/cbortoluzzi/wisent_project/annotation/WISENTM_Urano_F1_hap1/run_1/galba.gtf > /cluster/work/pausch/cbortoluzzi/wisent_project/annotation/WISENTM_Urano_F1_hap1/run_1/tsebra.log 2> /cluster/work/pausch/cbortoluzzi/wisent_project/annotation/WISENTM_Urano_F1_hap1/run_1/errors/tsebra.err
python3 /cluster/work/pausch/cbortoluzzi/softwares/Augustus-3.5.0/scripts/getAnnoFastaFromJoingenes.py -g /cluster/work/pausch/cbortoluzzi/wisent_project/repeat_annotation/Bison_bonasus/RepeatMasker/WISENTM_Urano_F1.hap1.rm.unmasked.fasta.masked -f /cluster/work/pausch/cbortoluzzi/wisent_project/annotation/WISENTM_Urano_F1_hap1/run_1/galba.gtf -o /cluster/work/pausch/cbortoluzzi/wisent_project/annotation/WISENTM_Urano_F1_hap1/run_1/galba 1>/cluster/work/pausch/cbortoluzzi/wisent_project/annotation/WISENTM_Urano_F1_hap1/run_1/getAnnoFastaFromJoingenes.tsebra.stdout 2> /cluster/work/pausch/cbortoluzzi/wisent_project/annotation/WISENTM_Urano_F1_hap1/run_1/errors/getAnnoFastaFromJoingenes.tsebra.stderr

echo -e "# IMPORTANT INFORMATION: the final output files of this GALBA run are galba.gtf, galba.aa, and galba.codingseq"
echo -e "This gene set is a result of running TSEBRA. In rare cases, the tsebra gene set may be too small due to a lack of evidence. In these cases, please compare to the augustus.hints.gtf gene set and use the one that is better"



