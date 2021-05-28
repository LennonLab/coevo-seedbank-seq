#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=02:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=SNCt-L1

##### load dependencies #####
module load samtools
module unload python
module load python/2.7.16

samtools mpileup -q10 -f /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/SNCt-L1-T14-pl/data/reference.fasta /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/SNCt-L1-T14-pl/data/reference.bam /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/SNCt-L1-T14-rV/data/reference.bam > /N/slate/danschw/coevo-seedbank-seq/data/CompSubPops/host/CompSubPop_merged/SNCt-L1.pileup

cat /N/slate/danschw/coevo-seedbank-seq/data/CompSubPops/host/CompSubPop_merged/SNCt-L1.pileup | python /N/slate/danschw/coevo-seedbank-seq/code/python/create_timecourse.py SNCt-L1 100 300 > /N/slate/danschw/coevo-seedbank-seq/data/CompSubPops/host/CompSubPop_merged/SNCt-L1_CompSubPop.txt
