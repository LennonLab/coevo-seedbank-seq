#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=2:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=WSO-L2-T14-rS

##### load dependencies #####
module load r
module load bowtie2

##### run breseq from local instance #####
/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin/breseq -j 8 -p --brief-html-output --polymorphism-reject-indel-homopolymer-length 0 --polymorphism-reject-surrounding-homopolymer-length 0 --polymorphism-score-cutoff 2 -o /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/unmatched/WSO-L2-T14-rS -r /N/slate/danschw/coevo-seedbank-seq/data/SPO1-ANC.gff3 /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq1/WSO-L2-T14-rS/data/ddup*R1*.fastq /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq1/WSO-L2-T14-rS/data/ddup*R2*.fastq > /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/unmatched_err/WSO-L2-T14-rS.err 2> /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/unmatched_log/WSO-L2-T14-rS.log