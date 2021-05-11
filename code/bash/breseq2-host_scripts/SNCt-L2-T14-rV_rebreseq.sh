#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=SNCt-L2-T14-rV

##### load dependencies #####
module load r
module load bowtie2

##### run breseq from local instance #####
/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin/breseq -j 8 -p --user-evidence-gd /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq_jc/merged/SNCt-L2.gd -o /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/SNCt-L2-T14-rV -r /N/slate/danschw/coevo-seedbank-seq/data/dSpoIIE-ANC2.gff3 /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L2-T14-rV_S51_R1_001.fastq.gz /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L2-T14-rV_S51_R2_001.fastq.gz > /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2_err/SNCt-L2-T14-rV.err 2> /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2_log/SNCt-L2-T14-rV.log
