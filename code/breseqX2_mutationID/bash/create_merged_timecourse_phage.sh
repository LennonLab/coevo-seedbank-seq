#!/bin/bash

# modified from 
# https://github.com/benjaminhgood/LTEE-metagenomic/blob/master/cluster_scripts/create_merged_timecourse.sh
# AND
# https://github.com/MURI2/Phylo_Evol_Timeseries/blob/master/bash/create_merged_timecourse.sh

module unload python
module load python/2.7.16 

##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
ODIR=${PARENT}/data/map-EVOL/phage
#path to python script that does the work
create_breseq_timecourse=${PARENT}/code/python/create_breseq_timecourse.py

mkdir -p ${PARENT}/data/timecourses/phage/timecourse_merged2

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
  gd_files="${ODIR}/breseq2/${pop}"-*"/output/evidence/evidence.gd"
    #all samples have the same reference, so only keep one
  ref="${ODIR}/breseq2/${pop}-T1/data/reference.fasta"

  merged_timecourse="${PARENT}/data/timecourses/phage/timecourse_merged2/${pop}_merged_timecourse.bz"
  if [ -f $merged_timecourse ]; then
    rm $merged_timecourse
  fi
 # write everything to a merged file
 cat "${PARENT}/data/timecourses/phage/timecourse_merged/${pop}_timecourse.txt" | python $create_breseq_timecourse $pop ${gd_files[@]} | bzip2 -c > $merged_timecourse
done



# write everything to a merged file
#cat data/timecourse_files/${population}_*timecourse.txt | python create_breseq_timecourse.py ${reference_file} ${population} ${gd_files} | bzip2 -c > ${population}_merged_timecourse.bz2
