#!/bin/bash

# modified from
# https://github.com/MURI2/Phylo_Evol_Timeseries/blob/master/bash/rebreseq.sh

# Script to identify mutations in sequenced host populations
# using breseq (polymorphism mode) to map reads to host ancestor



# path to local instance of breseq 0.35.5 
BRESEQ=/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin



##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
cd $PARENT/data

bash_breseq_scripts="${PARENT}/code/bash/unmatched-host_scripts"
data_breseq="${PARENT}/data/map-EVOL/host/unmatched"
#data_breseq_out="${PARENT}/data/map-EVOL/host/unmatched_out"
data_breseq_err="${PARENT}/data/map-EVOL/host/unmatched_err"
data_breseq_log="${PARENT}/data/map-EVOL/host/unmatched_log"

mkdir -p $bash_breseq_scripts
mkdir -p $data_breseq
#mkdir -p $data_breseq_out
mkdir -p $data_breseq_err
mkdir -p $data_breseq_log

# phage ancestor 
PREF=$PARENT/data/SPO1-ANC.gff3

# unmatched reads
dREADS=$PARENT/data/map-EVOL/host/breseq1




#####define population samples #####
declare -a TRTS=("WLO" "WSO" "SNO" )
declare -a REPS=("L1" "L2" "L3")
declare -a TIMES=("T1" "T4" "T14")
#declare -a TIMES=("T1" "T4" "T7" "T10" "T14")
declare -a EXTRACTS=("rV" "rS") # revived veg, revived spore

export -a SAMPLES=()

for TRT in "${TRTS[@]}"
do
	for REP in "${REPS[@]}"
	do
		for TIME in "${TIMES[@]}"
		do
			SAMPLES+=("${TRT}-${REP}-${TIME}-${EXTRACTS[0]}")
			SAMPLES+=("${TRT}-${REP}-${TIME}-${EXTRACTS[1]}")		
			if [ ${TRT} != "SNO" ]; then
				SAMPLES+=("${TRT}-${REP}-${TIME}-${EXTRACTS[2]}")
			fi
		done
	done
done


##### map samples #####
for SAMPLE in "${SAMPLES[@]}"
do
	# reads for current sample
	fR1=$dREADS/$SAMPLE/data/ddup*R1*.fastq
	fR2=$dREADS/$SAMPLE/data/ddup*R2*.fastq

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
    echo '#SBATCH --mem=50gb' >> $bash_out
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
