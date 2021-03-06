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

######  get host reference geneome B. subtilis delta6  #####

PARENT=/N/slate/danschw/coevo-seedbank-seq
cd $PARENT/data

#get host reference geneome (delta6)
wget  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/660/525/GCF_001660525.1_ASM166052v1/GCF_001660525.1_ASM166052v1_genomic.gbff.gz
#unzip ref genome
gunzip GCF*.gz

# for the spoIIE deletion strain we have manually edited the Delta6 genome
# the CDS of spoIIE (except for start and stop codons was replaced with the
# molecular scar left after cre/lox marker removal, as dexscribed by 
# Koo et al. 2017. Cell Systems 4:291


######  map sequenced ANC to ref #####

# Define paths

# host ancestor ref
wtHREF=$PARENT/data/GCF_001660525.1*
mutHREF=$PARENT/data/dSpoIIE_GCF_001660525.1*
# host ancestor reads
wtR1=$(find $PARENT/data/ddup-fastq/host/ -name "ddup*delta6-ANC*R1*.gz")
wtR2=$(find $PARENT/data/ddup-fastq/host/ -name "ddup*delta6-ANC*R2*.gz")

mutR1=$(find $PARENT/data/ddup-fastq/host/ -name "ddup*dSpoIIE-ANC*R1*.gz")
mutR2=$(find $PARENT/data/ddup-fastq/host/ -name "ddup*dSpoIIE-ANC*R2*.gz")

#make directory for results
wtODIR=$PARENT/data/map-ANC/host2delta6/delta6
mkdir -p ${wtODIR}
mutODIR=$PARENT/data/map-ANC/host2delta6/dSPOIIE
mkdir -p ${mutODIR}

# run breseq in clonal mode 
$BRESEQ/breseq -r $wtHREF -l 500 -j 12 -n "delta6_ANC" -o $wtODIR $wtR1 $wtR2

$BRESEQ/breseq -r $mutHREF -l 500 -j 12 -n "dSpoIIE_ANC" -o $mutODIR $mutR1 $mutR2

# -l limits the coverage


######  Generating mutated host reference sequences #####
# Generating a mutated reference sequences
# depends on breseq and dependencies loaded above

wtGD=$wtODIR/output/output.gd

$BRESEQ/gdtools APPLY -f GFF3 -o $PARENT/data/delta6-ANC.gff3 -r $wtHREF $wtGD

mutGD=$mutODIR/output/output.gd

$BRESEQ/gdtools APPLY -f GFF3 -o $PARENT/data/dSpoIIE-ANC.gff3 -r $mutHREF $mutGD
