#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=SNO-L1-T14-pl

##### load dependencies #####
module load r
module load bowtie2

##### run breseq from local instance #####
/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin/breseq -j 8 -p --user-evidence-gd /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq_jc/merged/SNO-L1.gd -o /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/SNO-L1-T14-pl -r /N/slate/danschw/coevo-seedbank-seq/data/dSpoIIE-ANC2.gff3 /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L1-T14-pl_S17_R1_001.fastq.gz /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L1-T14-pl_S17_R2_001.fastq.gz > /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2_err/SNO-L1-T14-pl.err 2> /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2_log/SNO-L1-T14-pl.log
