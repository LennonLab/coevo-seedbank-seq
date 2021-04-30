To run the analysis pipeline:

(1) place the read fastq.gz data in the data/input/fastq folder. use separate folders for host and phage:
Host files: coevo-seedbank-seq/data/input/fastq/host
Phage files: coevo-seedbank-seq/data/input/fastq/phage

(2) Remove PCR duplicates from sequence data using FastUniq v1.1 (https://sourceforge.net/projects/fastuniq/files/)
*The path to the Fastuniq installation should be specified t the head of code/R/rmv-pcr-dups.R*
Files are organized and prepared for analysis by Fastuniq by runing code/R/rmv-pcr-dups.R, from within the code folder:
$ module load r
$ Rscript code/R/rmv-pcr-dups.R host
$ Rscript code/R/rmv-pcr-dups.R phage

The R script writes a bash file with commands to run fastuniq in a slurm batch job, placed in coevo-seedbank-seq/code/bash.  From that folder run the commands:
$ sbatch deduplicate-host.sh
$ sbatch deduplicate-phage.sh

clean up after fastuniq:
$ cd coevo-seedbank-seq/data/ddup-fastq/
delete copied fastq files
$ find | grep -v "ddup" | grep "fastq" | xargs rm 
compress deduplicated files
$ find | grep "ddup" | xargs gzip 

(3) Construct ancestral reference genomes
Run bash script on slurm to map reads of ancestors to published reference
$ sbatch code/bash/map-anc-phage.sh
