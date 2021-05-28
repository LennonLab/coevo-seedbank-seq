#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=10gb,walltime=72:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

module unload python
module load python/3.6.1
python /N/slate/danschw/coevo-seedbank-seq/python/annotate_pvalues.py /N/slate/danschw/coevo-seedbank-seq/data/timecourses/phage/timecourse_depth/WLO-L1_depth_timecourse.bz /N/slate/danschw/coevo-seedbank-seq/data/timecourses/phage/timecourse_snp/WLO-L1_snp_timecourse.bz /N/slate/danschw/coevo-seedbank-seq/data/timecourses/phage/timecourse_indel/WLO-L1_indel_timecourse.bz /N/slate/danschw/coevo-seedbank-seq/data/timecourses/phage/timecourse_likelihood/WLO-L1_likelihood_timecourse.bz
