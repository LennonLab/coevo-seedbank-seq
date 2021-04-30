#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=23:59:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=Map-EvolPhage

# Script to identify mutations in sequences phage populations
# using breseq (polymorphism mode) to map reads to phage ancestor

##### load dependencies #####
module load r
module load bowtie2

# path to local instance of breseq 0.35.5 
BRESEQ=/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin



##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
cd $PARENT/data

# phage ancestor 
PREF=$PARENT/data/SPO1-ANC.gff3

# de-duplicated reads
READS=$PARENT/data/ddup-fastq/phage
#make directory for results (first breseq iteration)
mkdir -p $PARENT/data/map-EVOL/phage/breseq1
ODIR=$PARENT/data/map-EVOL/phage/breseq1




#####define population samples #####
declare -a TRTS=("WLO" "WSO" "SNO")
declare -a REPS=("L1" "L2" "L3")
declare -a TIMES=("T1" "T4" "T7" "T10" "T14")
declare -a SAMPLES=()

for TRT in "${TRTS[@]}"
do
	for REP in "${REPS[@]}"
	do
		for TIME in "${TIMES[@]}"
		do
			SAMPLES+=("${TRT}-${REP}-${TIME}")
		done
	done
done


##### map samples #####
for SAMPLE in "${SAMPLES[@]}"
	do
		# reads for current sample
		fR1=`find $READS -name "ddup*-$SAMPLE\_*R1*.gz"`
		fR2=`find $READS -name "ddup*-$SAMPLE\_*R2*.gz"`

		# make output directory
		mkdir $ODIR/$SAMPLE

		# map with breseq
		$BRESEQ/breseq -j 8 -p --brief-html-output --polymorphism-reject-indel-homopolymer-length 0 --polymorphism-reject-surrounding-homopolymer-length 0 --polymorphism-score-cutoff 2 -o $ODIR/$SAMPLE -r $PREF $fR1 $fR2 > $ODIR/$SAMPLE/breseq.err 2> $ODIR/$SAMPLE/breseq.log >> $ODIR/$SAMPLE/bash.out
done
