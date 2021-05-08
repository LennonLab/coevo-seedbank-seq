#!/bin/bash

module load trimmomatic/0.36

mkdir -p /N/slate/danschw/coevo-seedbank-seq/data/trim-test
cd /N/slate/danschw/coevo-seedbank-seq/data/trim-test

#nextera adapter file from the trimmomatic github site
# Trimmomatic could not read its own adapter files, even though they exist
# Having a local copy with a modifies name (smaal 'n') is the only way I got this thing to work
cp /N/slate/danschw/coevo-seedbank-seq/data/nexteraPE-PE.fa .

mkdir trimmed

#copy test files
cp /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/*delta6-ANC*.gz .
cp /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/*GSF2865-WLO-L1-T14-pl*.gz .

trimmomatic PE -phred33 -threads 1 \
/N/slate/danschw/coevo-seedbank-seq/data/trim-test/ddup-GSF2865-delta6-ANC-ANC-pl_S1_R1_001.fastq.gz \
/N/slate/danschw/coevo-seedbank-seq/data/trim-test/ddup-GSF2865-delta6-ANC-ANC-pl_S1_R2_001.fastq.gz \
/N/slate/danschw/coevo-seedbank-seq/data/trim-test/trimmed/trimmed-GSF2865-delta6-ANC-ANC-pl_S1_R1_001.fastq \
/N/slate/danschw/coevo-seedbank-seq/data/trim-test/trimmed/trimmed-GSF2865-delta6-ANC-ANC-pl_S1_R2_001.fastq \
ILLUMINACLIP:nexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> trim.paired_ANC.txt

trimmomatic PE -phred33 -threads 1 \
/N/slate/danschw/coevo-seedbank-seq/data/trim-test/ddup-GSF2865-WLO-L1-T14-pl_S5_R1_001.fastq.gz \
/N/slate/danschw/coevo-seedbank-seq/data/trim-test/ddup-GSF2865-WLO-L1-T14-pl_S5_R2_001.fastq.gz \
/N/slate/danschw/coevo-seedbank-seq/data/trim-test/trimmed/trimmed-ddup-GSF2865-WLO-L1-T14-pl_S5_R1_001.fastq \
/N/slate/danschw/coevo-seedbank-seq/data/trim-test/trimmed/trimmed-ddup-GSF2865-WLO-L1-T14-pl_S5_R2_001.fastq \
ILLUMINACLIP:nexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> trim.paired_WLO-L1-T14-pl.txt

