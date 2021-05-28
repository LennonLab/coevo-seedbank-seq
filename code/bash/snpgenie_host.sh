#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=host_snpGenie_test

##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
bash_scripts=${PARENT}/code/bash/snpgenie-host_scripts
ODIR=${PARENT}/data/SnpGenie/host
BRSQ=${PARENT}/data/map-EVOL/host/breseq2
SNPG=/N/u/danschw/Carbonate/my_tools/SNPGenie
GTF=/N/u/danschw/Carbonate/my_tools/gff_gtf/gffread

mkdir -p $bash_scripts
mkdir -p $ODIR
cd $ODIR

# iterate over all samples processed by breseq
declare -a SAMPLES=($(ls $BRSQ))

#declare -a SAMPLES="SNCt-L1-T14-pl"

for SAMPLE in "${SAMPLES[@]}";
do

	    bash_out="${bash_scripts}/${SAMPLE}_snpgenie.sh"
		    if [ -f $bash_out ]; then
		      rm -f $bash_out
		    fi
    
		echo '#!/bin/bash' >> $bash_out
		echo '#SBATCH --mail-user=danschw@iu.edu' >> $bash_out
		echo '#SBATCH --nodes=1' >> $bash_out
		echo '#SBATCH --ntasks-per-node=8' >> $bash_out
		echo '#SBATCH --time=3:59:00' >> $bash_out
		echo '#SBATCH --mem=50gb' >> $bash_out
		echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
		echo "#SBATCH --job-name=${SAMPLE}" >> $bash_out
		echo '' >> $bash_out
		echo "$GTF/gffread $BRSQ/$SAMPLE/data/reference.gff3 -T -o ${SAMPLE}_reference.gtf 2> ${SAMPLE}_gffread.log" >> $bash_out
		echo '# fix gtf file for SnpGenie' >> $bash_out
		echo "sed -i 's/transcript_id/gene_id/g' ${SAMPLE}_reference.gtf " >> $bash_out
		echo '' >> $bash_out
		echo '##### run snpgenie from local instance #####' >> $bash_out
		echo "perl $SNPG/snpgenie.pl  --minfreq=0.01 --snpreport=$BRSQ/$SAMPLE/data/output.vcf --vcfformat=2 --fastafile=$BRSQ/$SAMPLE/data/reference.fasta --gtffile=${SAMPLE}_reference.gtf --outdir $ODIR/$SAMPLE;" >> $bash_out
		echo '' >> $bash_out
		echo "rm ${SAMPLE}_reference.gtf ${SAMPLE}_gffread.log" >> $bash_out
	   	
		sbatch $bash_out

done


