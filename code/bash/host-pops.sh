#####define population samples #####
declare -a TRTS=("WLO" "WSO" "SNO" "WLCt" "WSCt" "SNCt")
declare -a REPS=("L1" "L2" "L3")
declare -a TIMES=("T14")
#declare -a TIMES=("T1" "T4" "T7" "T10" "T14")
declare -a EXTRAS=("pl" "rV" "rS")
export -a SAMPLES=()

# manually add the fonder populations
SAMPLES=("delta6-founder-T0-pl" "SpoIIE-founder-T0-pl")

for TRT in "${TRTS[@]}"
do
for REP in "${REPS[@]}"
do
for TIME in "${TIMES[@]}"
do
SAMPLES+=("${TRT}-${REP}-${TIME}-${EXTRAS[0]}")
SAMPLES+=("${TRT}-${REP}-${TIME}-${EXTRAS[1]}")		
if [ ${TRT} != "SNO" ]; then
SAMPLES+=("${TRT}-${REP}-${TIME}-${EXTRAS[2]}")
fi
done
done
done


#echo ${SAMPLES[@]}
 
for SAMPLE in "${SAMPLES[@]}"
do
#refernce for current sample
if [[ ${SAMPLE} == *"SN"* ]] || [[ ${SAMPLE} == *"SpoIIE"* ]] ; then
#echo "SSS-${SAMPLE}"
continue
else
#echo "WWW-${SAMPLE}"
SAMPLES=("${SAMPLES[@]/$SAMPLE}")
fi		
done

STR='GNU/Linux is an operating system'
SUB='Linux'
if [[ "$STR" == *"$SUB"* ]]; then
  echo "It's there."
fi















##### map samples #####
for SAMPLE in "${SAMPLES[@]}"
do
	# reads for current sample
	fR1=`find $dREADS -name "ddup*-$SAMPLE\_*R1*.gz"`
	fR2=`find $dREADS -name "ddup*-$SAMPLE\_*R2*.gz"`

	#refernce for current sample
	if [[ ${SAMPLE} == *"SN"* ]] || [[ ${SAMPLE} == *"SpoIIE"* ]] ; then
		HREF=$mutHREF
	else
		HREF=$wtHREF
	fi

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
    echo '#SBATCH --ntasks-per-node=8' >> $bash_out
    echo '#SBATCH --time=9:59:00' >> $bash_out
    echo '#SBATCH --mem=50gb' >> $bash_out
    echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
    echo "#SBATCH --job-name=${SAMPLE}" >> $bash_out
    echo '' >> $bash_out
    echo '##### load dependencies #####' >> $bash_out
    echo 'module load r' >> $bash_out
    echo 'module load bowtie2' >> $bash_out
    echo '' >> $bash_out
    echo '##### run breseq from local instance #####' >> $bash_out
    echo "$BRESEQ/breseq -j 8 -p --brief-html-output --polymorphism-reject-indel-homopolymer-length 0 --polymorphism-reject-surrounding-homopolymer-length 0 --polymorphism-score-cutoff 2 -o $OUT_breseq -r $HREF $fR1 $fR2 > $OUT_breseq_err 2> $OUT_breseq_log" >> $bash_out

   sbatch $bash_out

  else
    continue
  fi

done
