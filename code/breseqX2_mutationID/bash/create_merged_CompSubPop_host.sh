#!/bin/bash

# modified from 
# https://github.com/benjaminhgood/LTEE-metagenomic/blob/master/cluster_scripts/create_merged_timecourse.sh
# AND
# https://github.com/MURI2/Phylo_Evol_Timeseries/blob/master/bash/create_merged_timecourse.sh

module unload python
module load python/2.7.16 

##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
INDIR=${PARENT}/data/map-EVOL/host
ODIR=${PARENT}/data/CompSubPops/host
#path to python script that does the work
create_breseq_timecourse=${PARENT}/code/python/create_breseq_CompSubPop.py

mkdir -p ${ODIR}/CompSubPop_merged2

##### Define populations  #####
declare -a TRTS=("WLO" "WSO" "SNO" "WLCt" "WSCt" "SNCt")
declare -a REPS=("L1" "L2" "L3")
declare -a POPS=()

for TRT in "${TRTS[@]}"
do
	for REP in "${REPS[@]}"
	do
		POPS+=("${TRT}-${REP}")
	done
done

# To make the different extraction methods parsble as "time"
# we copy the gd files and change there names
old_gd_files=(${INDIR}/breseq2/*/output/evidence/evidence.gd)

for f in ${old_gd_files[@]}
do
	subpop="$(echo "$f" | sed s/.output.evidence.evidence.gd//g | sed s/.*breseq2.//g )"
	case "$(echo $subpop | sed s/.*-//g)" in
		"pl") subpop="$(echo $subpop | sed s/T14-pl/T100/g)";;
		"rS") subpop="$(echo $subpop | sed s/T14-rS/T200/g)";;
		"rV") subpop="$(echo $subpop | sed s/T14-rV/T300/g)";;
	esac
	cp $f ${ODIR}/CompSubPop_merged2/${subpop}.gd
done


for pop in "${POPS[@]}"
do
  gd_files="${ODIR}/CompSubPop_merged2/${pop}-*.gd"
    #all samples have the same reference, so only keep one
  ref="${INDIR}/breseq2/${pop}-T14-pl/data/reference.fasta"

  merged_comp="${ODIR}/CompSubPop_merged2/${pop}_merged_CompSubPop.bz"
  if [ -f $merged_comp ]; then
    rm $merged_comp
  fi
 # write everything to a merged file
 cat "${ODIR}/CompSubPop_merged/${pop}_CompSubPop.txt" | python $create_breseq_timecourse $pop ${gd_files[@]} | bzip2 -c > $merged_comp
  rm ${gd_files[@]}
done



# write everything to a merged file
#cat data/timecourse_files/${population}_*timecourse.txt | python create_breseq_timecourse.py ${reference_file} ${population} ${gd_files} | bzip2 -c > ${population}_merged_timecourse.bz2
