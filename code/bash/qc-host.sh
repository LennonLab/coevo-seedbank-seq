#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=QC-Hdedup

#### Load modules for fastQC ####

module load java
module unload perl/5.24.1
module load perl/5.30.1
module load fastqc/0.11.5 

##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
ODIR=${PARENT}/data/ddup-fastq/qc-host
mkdir -p ${ODIR}
#### fastQC ####

cd ${PARENT}/data/ddup-fastq/host

fastqc -o ${ODIR} ./*.fastq.gz

#### summarize with multiqc ####
module unload python
module load python/3.6.8
module load multiqc/python3.6/1.8

multiqc --interactive -o ${ODIR} ${ODIR} 
