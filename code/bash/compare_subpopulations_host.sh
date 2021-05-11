#!/bin/bash

# modified from 
# https://github.com/MURI2/Phylo_Evol_Timeseries/blob/master/bash/create_timecourse.sh

# we compare the 3 extraction sub-populations from each population
# "pl" = pellet, "rV" = revived veg, "rS" = revived spore
# taking atvantage of time course script, we encode these as integers

##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
#path to python script that does the work
create_timecourse=${PARENT}/code/python/create_timecourse.py
#path to results of second breseq
data_rebreseq="${PARENT}/data/map-EVOL/host/breseq2"
# make folders for current analysis
bash_CompSubPop_scripts="${PARENT}/code/bash/CompSubPop_scripts/host"
data_CompSubPop_out="${PARENT}/data/CompSubPops/host/CompSubPop_merged"

mkdir -p $data_CompSubPop_out
mkdir -p $bash_CompSubPop_scripts

##### Define populations  #####
#%declare -a TRTS=("WLO" "WSO" "SNO" "WLCt" "WSCt" "SNCt")
declare -a TRTS=("WLO" "WSO" "WLCt" "WSCt")
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
  bash_out="${bash_CompSubPop_scripts}/${pop}_CompSubPop.sh"
  if [ -f $bash_out ]; then
    rm $bash_out
  fi

  bam_files="${PARENT}/data/map-EVOL/host/breseq2/${pop}"-*"/data/reference.bam"

  declare -a TIMES=()
  #collect time in same order as bam files passed to pileup
    for pop_sample in ${bam_files[@]}
	do
    # extract subpopulation from bam file path
    subpop="$(echo "$pop_sample" | sed s/.data.reference.bam//g | sed s/.*-//g )"
		#convert to integer
		case $subpop in
			"pl") time=100;;
			"rS") time=200;;
			"rV") time=300;;
		esac

    
		TIMES+=("${time}")
  done

 
  #all samples have the same reference, so only keep one
  ref="${PARENT}/data/map-EVOL/host/breseq2/${pop}-T14-pl/data/reference.fasta"

  out="${data_CompSubPop_out}/${pop}.pileup"
  if [ -f $out ]; then
    rm $out
  fi

  out_CompSubPop="${data_CompSubPop_out}/${pop}_CompSubPop.txt"
  if [ -f $out_CompSubPop ]; then
    rm $out_CompSubPop
  fi

  echo '#!/bin/bash' >> $bash_out
  echo '#SBATCH --mail-user=danschw@iu.edu' >> $bash_out
  echo '#SBATCH --nodes=1' >> $bash_out
  echo '#SBATCH --ntasks-per-node=8' >> $bash_out
  echo '#SBATCH --time=02:00:00' >> $bash_out
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
  echo "cat ${out} | python ${create_timecourse} ${pop} ${TIMES[@]} > ${out_CompSubPop}" >> $bash_out

	
  echo "${pop}"
  echo "${TIMES[@]}"

  sbatch $bash_out
done



#for file in .* *;
#  do
#    echo $file
#    bzcat $file | head -10;
#done
