#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=MapAncHost

# load dependencies
module load r
module load bowtie2

# path to local instance of breseq 0.35.5 
BRESEQ=/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin

######  get host reference geneome B. subtilis 168  #####
# using this rather than Delta6 for better annotations
PARENT=/N/slate/danschw/coevo-seedbank-seq
cd $PARENT/data

#get host reference geneome (SPO1)
wget  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.gbff.gz

#unzip ref genome
gunzip GCF*.gz



######  map sequenced ANC to ref #####

# Define paths

# host ancestor ref
HREF=$PARENT/data/GCF_000009045.1*
# host ancestor reads
wtR1=$(find $PARENT/data/ddup-fastq/host/ -name "ddup*delta6-ANC*R1*.gz")
wtR2=$(find $PARENT/data/ddup-fastq/host/ -name "ddup*delta6-ANC*R2*.gz")

mutR1=$(find $PARENT/data/ddup-fastq/host/ -name "ddup*dSpoIIE-ANC*R1*.gz")
mutR2=$(find $PARENT/data/ddup-fastq/host/ -name "ddup*dSpoIIE-ANC*R2*.gz")

#make directory for results
wtODIR=$PARENT/data/map-ANC/host/delta6
mkdir -p ${wtODIR}
mutODIR=$PARENT/data/map-ANC/host/dSPOIIE
mkdir -p ${mutODIR}

# run breseq in clonal mode 
$BRESEQ/breseq -r $HREF -l 500 -j 12 -n "delta6_ANC" -o $wtODIR $wtR1 $wtR2

$BRESEQ/breseq -r $HREF -l 500 -j 12 -n "dSpoIIE_ANC" -o $mutODIR $mutR1 $mutR2

# -l limits the coverage


######  Generating mutated host reference sequences #####
# Generating a mutated reference sequences
# depends on breseq and dependencies loaded above

wtGD=$wtODIR/output/output.gd

$BRESEQ/gdtools APPLY -f GFF3 -o $PARENT/data/delta6-ANC.gff3 -r $HREF $wtGD

mutGD=$mutODIR/output/output.gd

$BRESEQ/gdtools APPLY -f GFF3 -o $PARENT/data/dSpoIIE-ANC.gff3 -r $HREF $mutGD
