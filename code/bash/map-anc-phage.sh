#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=Map-AncPhage

# load dependencies
module load r
module load bowtie2

# path to local instance of breseq 0.35.5 
BRESEQ=/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin

######  get phage reference geneome (SPO1)  #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
cd $PARENT/data

#get phage reference geneome (SPO1)
wget  wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Bacillus_virus_SPO1/latest_assembly_versions/GCF_000881675.1_ViralProj32379/GCF_000881675.1_ViralProj32379_genomic.gbff.gz

#unzip ref genome
gunzip GCF*.gz



######  map sequenced ANC to ref #####

# Define paths

# phage ancestor ref
PREF=$PARENT/data/GCF_000881675.1*
# phage ancestor reads
R1=$(find $PARENT/data/ddup-fastq/phage/ -name "ddup*ANC*R1*.gz")
R2=$(find $PARENT/data/ddup-fastq/phage/ -name "ddup*ANC*R2*.gz")

#make directory for results
mkdir -p $PARENT/data/map-ANC/phage
ODIR=$PARENT/data/map-ANC/phage

# run breseq in clonal mode 
$BRESEQ/breseq -r $PREF -l 500 -j 8 -n "SPO1_ANC" -o $ODIR $R1 $R2

# -l limits the coverage


######  Generating a mutated phage reference sequence #####
# Generating a mutated phage reference sequence
# depends on breseq and dependencies loaded above

GD=$PARENT/data/map-ANC/phage/output/output.gd

$BRESEQ/gdtools APPLY -f GFF3 -o $PARENT/data/SPO1-ANC.gff3 -r $PREF $GD
