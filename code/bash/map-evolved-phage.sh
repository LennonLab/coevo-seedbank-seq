#!/bin/bash

# Script to identify mutations in sequences phage populations
# using breseq (polymorphism mode) to map reads to phage ancestor

##### load dependencies #####

# path to local instance of breseq 0.36.1 
BRESEQ=/N/u/danschw/Carbonate/my_tools/breseq-0.36.1-Linux-x86_64/bin



##### Define paths #####
PARENT=/N/slate/danschw/Github/coevo-seedbank-seq
cd $PARENT/data

bash_breseq_scripts="${PARENT}/code/bash/breseq1-phage_scripts"
data_breseq="${PARENT}/data/map-EVOL/phage/breseq1"
data_breseq_err="${PARENT}/data/map-EVOL/phage/breseq1_err"
data_breseq_log="${PARENT}/data/map-EVOL/phage/breseq1_log"

mkdir -p $bash_breseq_scripts
mkdir -p $data_breseq
#mkdir -p $data_breseq_out
mkdir -p $data_breseq_err
mkdir -p $data_breseq_log

# phage ancestor 
PREF=$PARENT/data/SPO1-ANC.gff3

# de-duplicated reads
dREADS=$PARENT/data/ddup-fastq/phage
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
	fR1=`find $dREADS -name "ddup*-$SAMPLE\_*R1*.gz"`
	fR2=`find $dREADS -name "ddup*-$SAMPLE\_*R2*.gz"`


  files_test=( ${fR1} ${fR2} )
  if [ -e "${files_test[0]}" ] && [ -e "${files_test[1]}" ]; then

    bash_out="${bash_breseq_scripts}/${SAMPLE}_breseq.sh"
    if [ -f $bash_out ]; then
      rm -f $bash_out
    fi

    OUT_breseq="${data_breseq}/${SAMPLE}"
    mkdir -p $OUT_breseq

    OUT_breseq_err="${data_breseq_err}/${SAMPLE}.err"
    if [ -f $OUT_breseq_err ]; then
      rm -f $OUT_breseq_err
    fi

    OUT_breseq_log="${data_breseq_log}/${SAMPLE}.log"
    if [ -f $OUT_breseq_log ]; then
      rm -f $OUT_breseq_log
    fi

    echo '#!/bin/bash' >> $bash_out
    echo '#SBATCH --mail-user=danschw@iu.edu' >> $bash_out
    echo '#SBATCH --nodes=1' >> $bash_out
    echo '#SBATCH --ntasks-per-node=1' >> $bash_out
    echo '#SBATCH --cpus-per-task=8' >> $bash_out
    echo '#SBATCH --time=2:59:00' >> $bash_out
    echo '#SBATCH --mem=8gb' >> $bash_out
    echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
    echo "#SBATCH --job-name=${SAMPLE}" >> $bash_out
    echo '' >> $bash_out
    echo '##### load dependencies #####' >> $bash_out
    echo 'module load r' >> $bash_out
    echo 'module load bowtie2' >> $bash_out
    echo '' >> $bash_out
    echo '##### run breseq from local instance #####' >> $bash_out
    echo "$BRESEQ/breseq -j 8 -p --brief-html-output --polymorphism-reject-indel-homopolymer-length 0 --polymorphism-reject-surrounding-homopolymer-length 0 --polymorphism-score-cutoff 2 -o $OUT_breseq -r $PREF $fR1 $fR2 > $OUT_breseq_err 2> $OUT_breseq_log" >> $bash_out

   sbatch $bash_out

  else
    continue
  fi

done

