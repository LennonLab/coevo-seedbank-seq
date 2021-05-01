#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=WLO-L1-T4

##### load dependencies #####
module load r
module load bowtie2

##### re-run breseq from local instance #####
/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin/breseq -j 8 -p --user-evidence-gd /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/phage/breseq_jc/merged/WLO-L1.gd -o /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/phage/breseq2/WLO-L1-T4 -r /N/slate/danschw/coevo-seedbank-seq/data/SPO1-ANC.gff3 /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-10-WLO-L1-T4_S48_R1_001.fastq.gz /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-10-WLO-L1-T4_S48_R2_001.fastq.gz > /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/phage/breseq2_out/WLO-L1-T4.out 2> /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/phage/breseq2_err/WLO-L1-T4.err
