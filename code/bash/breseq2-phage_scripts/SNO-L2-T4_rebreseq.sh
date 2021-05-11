#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=SNO-L2-T4

##### load dependencies #####
module load r
module load bowtie2

##### re-run breseq from local instance #####
/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin/breseq -j 8 -p --user-evidence-gd /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/phage/breseq_jc/merged/SNO-L2.gd -o /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/phage/breseq2/SNO-L2-T4 -r /N/slate/danschw/coevo-seedbank-seq/data/SPO1-ANC.gff3 /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-17-SNO-L2-T4_S30_R1_001.fastq.gz /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-17-SNO-L2-T4_S30_R2_001.fastq.gz > /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/phage/breseq2_out/SNO-L2-T4.out 2> /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/phage/breseq2_err/SNO-L2-T4.err
