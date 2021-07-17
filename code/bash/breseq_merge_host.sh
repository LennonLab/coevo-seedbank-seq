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
ODIR=${PARENT}/data/map-EVOL/host

mkdir -p ${ODIR}/breseq_merge/
mkdir -p ${ODIR}/breseq_merge/all
mkdir -p ${ODIR}/breseq_merge/merged



#####define population samples #####
declare -a TRTS=("WLO" "WSO" "SNO" "WLCt" "WSCt" "SNCt")
declare -a REPS=("L1" "L2" "L3")
declare -a TIMES=("T1" "T4" "T14")
#declare -a TIMES=("T1" "T4" "T7" "T10" "T14")
declare -a SUBPOPS=("rV" "rS") # revived veg, revived spore
declare -a SAMPLES=()

# add founder populations to all
#wt_evidence_founder="${ODIR}/breseq1/delta6-founder-T0-pl/output/evidence/evidence.gd"
#mut_evidence_founder="${ODIR}/breseq1/dSpoIIE-founder-T0-pl/output/evidence/evidence.gd"

##### merge evidence files of all samples per population #####

for TRT in "${TRTS[@]}"
do
  for REP in "${REPS[@]}"
  do
      
    POP=${TRT}-${REP}
    declare -a gd_files=()
    
    for SUBPOP in "${SUBPOPS[@]}"
    do
      for TIME in "${TIMES[@]}"
      do
        
        current_gd="${ODIR}/breseq1/${TRT}-${REP}-${TIME}-${SUBPOP}/output/evidence/evidence.gd"
        if [ -f $current_gd ] ; then
					gd_files=(${gd_files[@]} $current_gd)
				fi

      done
    done
  
      # make new file
			merged_output="${ODIR}/breseq_merge/merged/${POP}.gd"
								if [ -f $merged_output ] ; then
									rm $merged_output
								fi
			${BRESEQ}/gdtools UNION -o $merged_output -e ${gd_files[@]}
			
  done
done

