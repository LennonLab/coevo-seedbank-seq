#!/bin/bash

# modified from
# https://github.com/MURI2/Bacillus_Evol_Timeseries/blob/master/bin/breseq_jc.sh

##### load dependencies #####
module load r
module load bowtie2

# path to local instance of breseq 0.35.5 
BRESEQ=/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin

##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
ODIR=${PARENT}/data/map-EVOL/phage/

mkdir -p ${ODIR}/breseq_jc/
mkdir -p ${ODIR}/breseq_jc/all
mkdir -p ${ODIR}/breseq_jc/merged



#####define population samples #####
declare -a TRTS=("WLO" "WSO" "SNO")
declare -a REPS=("L1" "L2" "L3")
declare -a TIMES=("T1" "T4" "T7" "T10" "T14")
declare -a SAMPLES=()



##### merge evidence files of all times per population #####

for TRT in "${TRTS[@]}"
do
	for REP in "${REPS[@]}"
	do
	      POP="${TRT}-${REP}"
	      for TIME in "${TIMES[@]}"
	      do
		evidence_time="${ODIR}/breseq1/${POP}-${TIME}/output/evidence/evidence.gd"
		if [ -f $evidence_time ] ; then
		  junction_output="${ODIR}/breseq_jc/all/${POP}-${TIME}.gd"
		  cat $evidence_time | grep 'JC\|#' > $junction_output
		fi
	      done
	      junction_merged_output="${ODIR}/breseq_jc/merged/${POP}.gd"
	      if [ -f $junction_merged_output ] ; then
		rm $junction_merged_output
	      fi
	      evidence="${ODIR}/breseq_jc/all/${POP}"*".gd"
	      ${BRESEQ}/gdtools UNION -o $junction_merged_output -e $evidence
	done
done
