#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=02:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=WLCt-L1

##### load dependencies #####
module load samtools
module unload python
module load python/2.7.16

samtools mpileup -q10 -f /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/WLCt-L1-T14-pl/data/reference.fasta /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/WLCt-L1-T14-pl/data/reference.bam /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/WLCt-L1-T14-rS/data/reference.bam /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/WLCt-L1-T14-rV/data/reference.bam > /N/slate/danschw/coevo-seedbank-seq/data/CompSubPops/host/CompSubPop_merged/WLCt-L1.pileup

cat /N/slate/danschw/coevo-seedbank-seq/data/CompSubPops/host/CompSubPop_merged/WLCt-L1.pileup | python /N/slate/danschw/coevo-seedbank-seq/code/python/create_timecourse.py WLCt-L1 100 200 300 > /N/slate/danschw/coevo-seedbank-seq/data/CompSubPops/host/CompSubPop_merged/WLCt-L1_CompSubPop.txt
