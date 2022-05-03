#!/bin/bash

# modified from
# https://github.com/MURI2/Bacillus_Evol_Timeseries/blob/master/bin/breseq_merge.sh

##### load dependencies #####
module load r
module load bowtie2

# path to local instance of breseq 0.35.5 
BRESEQ=/N/u/danschw/Carbonate/my_tools/breseq-0.35.5-Linux-x86_64/bin

##### Define paths #####
PARENT=/N/slate/danschw/coevo-seedbank-seq
ODIR=${PARENT}/data/map-EVOL/phage/

mkdir -p ${ODIR}/breseq_merge/




#####define population samples #####
declare -a TRTS=("WLO" "WSO" "SNO")
declare -a REPS=("L1" "L2" "L3")
declare -a TIMES=("T1" "T4" "T7" "T10" "T14")
declare -a SAMPLES=()



##### merge evidence files of all samples per host #####
declare -a phage_gd_files=()

for TRT in "${TRTS[@]}"
do
  for REP in "${REPS[@]}"
  do
    POP="${TRT}-${REP}"
    for TIME in "${TIMES[@]}"
      do
        current_gd="${ODIR}/breseq1/${TRT}-${REP}-${TIME}/output/evidence/evidence.gd"
        if [ -f $current_gd ] ; then
          phage_gd_files=(${phage_gd_files[@]} $current_gd)
        fi
    done
  done
done

merged_output="${ODIR}/breseq_merge/phage-pops-merged.gd"
					if [ -f $merged_output ] ; then
						rm $merged_output
					fi
${BRESEQ}/gdtools UNION -o $merged_output -e ${phage_gd_files[@]}