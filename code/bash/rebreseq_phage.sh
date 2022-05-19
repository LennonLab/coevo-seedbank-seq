#!/bin/bash

# modified from
# https://github.com/MURI2/Phylo_Evol_Timeseries/blob/master/bash/rebreseq.sh

##### Define paths #####
PARENT=/N/slate/danschw/Github/coevo-seedbank-seq
# de-duplicated reads
dREADS=${PARENT}/data/ddup-fastq/phage

bash_rebreseq_scripts="${PARENT}/code/bash/breseq2-phage_scripts"
data_rebreseq="${PARENT}/data/map-EVOL/phage/breseq2"
data_rebreseq_err="${PARENT}/data/map-EVOL/phage/breseq2_err"
data_rebreseq_log="${PARENT}/data/map-EVOL/phage/breseq2_log"


mkdir -p $bash_rebreseq_scripts
mkdir -p $data_rebreseq
mkdir -p $data_rebreseq_log
mkdir -p $data_rebreseq_err

# phage ancestor 
PREF=${PARENT}/data/SPO1-ANC.gff3

# path to local instance of breseq 0.36.1 
BRESEQ=/N/u/danschw/Carbonate/my_tools/breseq-0.36.1-Linux-x86_64/bin

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



for SAMPLE in "${SAMPLES[@]}"
do
	# reads for current sample
	fR1=`find $dREADS -name "ddup*-$SAMPLE\_*R1*.gz"`
	fR2=`find $dREADS -name "ddup*-$SAMPLE\_*R2*.gz"`

  files_test=( ${fR1} ${fR2} )
  if [ -e "${files_test[0]}" ] && [ -e "${files_test[1]}" ]; then
    pop="$(echo "$SAMPLE" | cut -d "-" -f1-2)"
    pop_gd="${PARENT}/data/map-EVOL/phage/breseq_merge/phage-pops-merged.gd"

    bash_out="${bash_rebreseq_scripts}/${SAMPLE}_rebreseq.sh"
    if [ -f $bash_out ]; then
      rm -f $bash_out
    fi
    OUT_rebreseq="${data_rebreseq}/${SAMPLE}"
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
    echo '#SBATCH --time=1:59:00' >> $bash_out
    echo '#SBATCH --mem=4gb' >> $bash_out
    echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
    echo "#SBATCH --job-name=p${SAMPLE}" >> $bash_out
    echo '' >> $bash_out
    echo '##### load dependencies #####' >> $bash_out
    echo 'module load r' >> $bash_out
    echo 'module load bowtie2' >> $bash_out
    echo '' >> $bash_out
    echo '##### re-run breseq from local instance #####' >> $bash_out
    echo "${BRESEQ}/breseq -j 8 -p --user-evidence-gd ${pop_gd} -o ${OUT_rebreseq} -r ${PREF} ${fR1} ${fR2} > ${OUT_rebreseq_err} 2> ${OUT_rebreseq_log}" >> $bash_out

    sbatch $bash_out

  else
    continue
  fi

done
