#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=ddup-WLCt-L2-T14-rS

#copy and extract paired read data
zcat /N/slate/danschw/coevo-seedbank-seq/data/input/fastq/host/GSF2865-WLCt-L2-T14-rS_S33_R1_001.fastq.gz > /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/GSF2865-WLCt-L2-T14-rS_S33_R1_001.fastq
zcat /N/slate/danschw/coevo-seedbank-seq/data/input/fastq/host/GSF2865-WLCt-L2-T14-rS_S33_R2_001.fastq.gz > /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/GSF2865-WLCt-L2-T14-rS_S33_R2_001.fastq

# make input list for FastUniq
touch /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/input_list_WLCt-L2-T14-rS.txt
echo GSF2865-WLCt-L2-T14-rS_S33_R1_001.fastq >> /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/input_list_WLCt-L2-T14-rS.txt
echo GSF2865-WLCt-L2-T14-rS_S33_R2_001.fastq >> /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/input_list_WLCt-L2-T14-rS.txt

#de-duplicate with FastUniq (local)
cd  /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host

 /N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/input_list_WLCt-L2-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L2-T14-rS_S33_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L2-T14-rS_S33_R2_001.fastq -c 1

 # delete copied inputs
rm  /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/GSF2865-WLCt-L2-T14-rS_S33_R1_001.fastq
rm  /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/GSF2865-WLCt-L2-T14-rS_S33_R2_001.fastq

 # delete input list
rm  /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/input_list_WLCt-L2-T14-rS.txt

 # compress outputs
gzip  /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L2-T14-rS_S33_R1_001.fastq
gzip  /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L2-T14-rS_S33_R2_001.fastq
