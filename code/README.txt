To run the analysis pipeline:

(1) place the read fastq.gz data in the data/input/fastq folder. use separate folders for host and phage:
Host files: coevo-seedbank-seq/data/input/fastq/host
Phage files: coevo-seedbank-seq/data/input/fastq/phage

(2) Remove PCR duplicates from sequence data using FastUniq v1.1 (https://sourceforge.net/projects/fastuniq/files/)
*The path to the Fastuniq installation should be specified at the head of code/R/rmv-pcr-dups.R*
The following code writes and submits slurm batch jobs for deuplicating.
Needs to be run separetly for phage and host
$ module load r
$ Rscript code/R/rmv-pcr-dups.R host
$ Rscript code/R/rmv-pcr-dups.R phage

Run FastQC and summarize with multiQC
$ sbatch code/bash/qc-host.sh
$ sbatch code/bash/qc-phage.sh

(tested need for trimming host to remove adapter sequences, this does not seem to have any effect. Skipping. to run test code: bash code/bash/trim_adapters-test.sh)

(3) Construct ancestral reference genomes
Run bash script on slurm to map reads of ancestors to published reference
$ sbatch code/bash/map-anc-phage.sh
$ sbatch code/bash/map-anc-hosts.sh

(4) Lists data per sample
List of al the samples to be used in downstream analysis
$ Rscript code/R/sample-list.R
Writes csv files in data folder (samples_evolved.csv; samples_ancestors.csv)

