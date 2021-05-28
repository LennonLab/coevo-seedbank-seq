#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=3:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=WLO-L1-T14-rV

/N/u/danschw/Carbonate/my_tools/gff_gtf/gffread/gffread /N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/WLO-L1-T14-rV/data/reference.gff3 -T -o reference.gtf 2> gffread.log
# fix gtf file for SnpGenie
sed -i 's/transcript_id/gene_id/g' reference.gtf 

##### run snpgenie from local instance #####
perl /N/u/danschw/Carbonate/my_tools/SNPGenie/snpgenie.pl  --minfreq=0.01 --snpreport=/N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/WLO-L1-T14-rV/data/output.vcf --vcfformat=2 --fastafile=/N/slate/danschw/coevo-seedbank-seq/data/map-EVOL/host/breseq2/WLO-L1-T14-rV/data/reference.fasta --gtffile=reference.gtf --outdir /N/slate/danschw/coevo-seedbank-seq/data/SnpGenie/host/WLO-L1-T14-rV;

rm reference.gtf gffread.log
