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

bash_rebreseq_scripts="${PARENT}/code/bash/breseq2-host_scripts"
data_breseq="${PARENT}/data/map-EVOL/host/breseq2"
#data_rebreseq_out="${PARENT}/data/map-EVOL/host/breseq2_out"
data_rebreseq_err="${PARENT}/data/map-EVOL/host/breseq2_err"
data_rebreseq_log="${PARENT}/data/map-EVOL/host/breseq2_log"

mkdir -p $bash_rebreseq_scripts
mkdir -p $data_breseq
#mkdir -p $data_rebreseq_out
mkdir -p $data_rebreseq_err
mkdir -p $data_rebreseq_log

# host ancestor 
wtHREF=$PARENT/data/delta6-ANC2.gff3
mutHREF=$PARENT/data/dSpoIIE-ANC2.gff3

# de-duplicated reads
dREADS=$PARENT/data/ddup-fastq/host
#make directory for results (first breseq iteration)
ODIR=$PARENT/data/map-EVOL/host/breseq2
mkdir -p $ODIR



#####define population samples #####
declare -a TRTS=("WLO" "WSO" "SNO" "WLCt" "WSCt" "SNCt")
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
	fR1=`find $dREADS -name "ddup*-$SAMPLE\_*R1*.gz"`
	fR2=`find $dREADS -name "ddup*-$SAMPLE\_*R2*.gz"`

	#refernce for current sample
	if [[ ${SAMPLE} == *"SN"* ]] ; then
		HREF=$mutHREF
	else
		HREF=$wtHREF
	fi

	files_test=( ${fR1} ${fR2} )
	if [ -e "${files_test[0]}" ] && [ -e "${files_test[1]}" ]; then

    	pop="$(echo "$SAMPLE" | cut -d "-" -f1-2)"
	    pop_gd="${PARENT}/data/map-EVOL/host/breseq_merge/merged/${pop}.gd"

    bash_out="${bash_rebreseq_scripts}/${SAMPLE}_rebreseq.sh"
    if [ -f $bash_out ]; then
      rm -f $bash_out
    fi

    OUT_rebreseq="${data_breseq}/${SAMPLE}"
    mkdir -p $OUT_rebreseq

    OUT_rebreseq_err="${data_rebreseq_err}/${SAMPLE}.err"
    if [ -f $OUT_rebreseq_err ]; then
      rm -f $OUT_rebreseq_err
    fi

    OUT_rebreseq_log="${data_rebreseq_log}/${SAMPLE}.log"
    if [ -f $OUT_rebreseq_log ]; then
      rm -f $OUT_rebreseq_log
    fi

    echo '#!/bin/bash' >> $bash_out
    echo '#SBATCH --mail-user=danschw@iu.edu' >> $bash_out
    echo '#SBATCH --nodes=1' >> $bash_out
    echo '#SBATCH --ntasks-per-node=1' >> $bash_out
    echo '#SBATCH --cpus-per-task=8' >> $bash_out
    echo '#SBATCH --time=2:59:00' >> $bash_out
    echo '#SBATCH --mem=20gb' >> $bash_out
    echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
    echo "#SBATCH --job-name=${SAMPLE}" >> $bash_out
    echo '' >> $bash_out
    echo '##### load dependencies #####' >> $bash_out
    echo 'module load r' >> $bash_out
    echo 'module load bowtie2' >> $bash_out
    echo '' >> $bash_out
    echo '##### run breseq from local instance #####' >> $bash_out
    echo "${BRESEQ}/breseq -j 8 -p --user-evidence-gd ${pop_gd} -o ${OUT_rebreseq} -r ${HREF} ${fR1} ${fR2} > ${OUT_rebreseq_err} 2> ${OUT_rebreseq_log}" >> $bash_out
   	sbatch $bash_out

  else
    continue
  fi

done
