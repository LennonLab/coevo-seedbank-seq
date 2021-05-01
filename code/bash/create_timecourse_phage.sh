#!/bin/bash

# modified from 
# https://github.com/MURI2/Phylo_Evol_Timeseries/blob/master/bash/create_timecourse.sh

##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
#path to python script that does the work
create_timecourse=${PARENT}/code/python/create_timecourse.py
#path to results of second breseq
data_rebreseq="${PARENT}/data/map-EVOL/phage/breseq2"
# make folders for current analysis
bash_timecourse_scripts="${PARENT}/code/bash/timecourse_scripts"
data_timecourse_out="${PARENT}/data/timecourses/phage/timecourse_merged"

mkdir -p $data_timecourse_out
mkdir -p $bash_timecourse_scripts

##### Define populations  #####
declare -a TRTS=("WLO" "WSO" "SNO")
declare -a REPS=("L1" "L2" "L3")
declare -a POPS=()

for TRT in "${TRTS[@]}"
do
	for REP in "${REPS[@]}"
	do
		POPS+=("${TRT}-${REP}")
	done
done



for pop in "${POPS[@]}"
do
  bash_out="${bash_timecourse_scripts}/${pop}_timecourse.sh"
  if [ -f $bash_out ]; then
    rm $bash_out
  fi

  declare -a times=()
  for pop_sample in `find ${data_rebreseq} -type d -name  *${pop}*`
  do
    #time="$(echo "$pop_sample" | cut -d "/" -f10-10 | cut -d "_" -f2-2)"
    time="$(echo $pop_sample | sed 's|.*-T||')"
    times+=("${time}")
  done
	# sort time array (not sure if required)
	TIMES="$(echo ${times[@]} | xargs -n1| sort -n | cut -d " " -f1-100 |xargs)"
	TIMES=($TIMES)


  bam_files="${PARENT}/data/map-EVOL/phage/breseq2/${pop}"-*"/data/reference.bam"
  
  #all samples have the same reference, so only keep one
  ref="${PARENT}/data/map-EVOL/phage/breseq2/${pop}-T1/data/reference.fasta"

  out="${data_timecourse_out}/${pop}.pileup"
  if [ -f $out ]; then
    rm $out
  fi

  out_timecourse="${data_timecourse_out}/${pop}_timecourse.txt"
  if [ -f $out_timecourse ]; then
    rm $out_timecourse
  fi

  echo '#!/bin/bash' >> $bash_out
  echo '#SBATCH --mail-user=danschw@iu.edu' >> $bash_out
  echo '#SBATCH --nodes=1' >> $bash_out
  echo '#SBATCH --ntasks-per-node=8' >> $bash_out
  echo '#SBATCH --time=9:59:00' >> $bash_out
  echo '#SBATCH --mem=50gb' >> $bash_out
  echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
  echo "#SBATCH --job-name=${pop}" >> $bash_out
  echo '' >> $bash_out
  echo '##### load dependencies #####' >> $bash_out
  echo 'module load samtools' >> $bash_out
  echo 'module unload python' >> $bash_out
  echo 'module load python/2.7.16' >> $bash_out
  echo '' >> $bash_out
  echo "samtools mpileup -q10 -f" ${ref} ${bam_files} ">" ${out} >> $bash_out
  echo '' >> $bash_out
  echo "cat ${out} | python ${create_timecourse} ${pop} ${TIMES[@]} > ${out_timecourse}" >> $bash_out

	
  echo "${pop}"
  echo "${TIMES[@]}"

  sbatch $bash_out
done



#for file in .* *;
#  do
#    echo $file
#    bzcat $file | head -10;
#done
