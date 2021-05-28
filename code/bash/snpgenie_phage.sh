#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=phage_snpGenie

##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq

ODIR=${PARENT}/data/SnpGenie/phage
BRSQ=${PARENT}/data/map-EVOL/phage/breseq2
SNPG=/N/u/danschw/Carbonate/my_tools/SNPGenie
GTF=/N/u/danschw/Carbonate/my_tools/gff_gtf/gffread

mkdir -p $ODIR
cd $ODIR

# iterate over all samples processed by breseq
declare -a SAMPLES=($(ls $BRSQ))

for SAMPLE in "${SAMPLES[@]}";
	do


		$GTF/gffread $BRSQ/$SAMPLE/data/reference.gff3 -T -o reference.gtf 2> gffread.log
		# fix gtf file for SnpGenie
		sed -i 's/transcript_id/gene_id/g' reference.gtf 

		perl $SNPG/snpgenie.pl  --minfreq=0.01 \
		--snpreport=$BRSQ/$SAMPLE/data/output.vcf --vcfformat=2 \
		--fastafile=$BRSQ/$SAMPLE/data/reference.fasta --gtffile=reference.gtf \
		--outdir $ODIR/$SAMPLE;

		rm reference.gtf gffread.log
done


