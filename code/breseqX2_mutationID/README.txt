# Old code for Ben Good type analyses
#===================================================
Cleaning up I stashed aside all the analyses code in this folder
(4) Initial breseq round. Map reads of evolved samples to ANC
$sbatch map-evolved-phage.sh
$bash map-evolved-hosts.sh
(5) Compile merged list of candidate junctions per population. 
$ bash breseq_jc_phage.sh
$ bash breseq_jc_host.sh

(6) Second round of breseq. Map reads of evolved samples to ANC 
The following bash script writes and submits sbatch jobs for all samples
$ bash rebreseq_phage.sh
$ bash rebreseq_host.sh

>>> mapping unmatched reads to phage: code/bash/map-umatched-hosts.sh
>>> anlysis of breseq2 coverage: code/R/breseq-json-coverage.R


(7) pileup and merge time courses
$ bash create_timecourse_phage.sh

$ compare_subpopulations_host.sh


(8) Create junction timecourses, and final merged timecourse file per
population
$ bash create_merged_timecourse_phage.sh

(9) make annotate_pvalues
based on BG's response to WRS (https://github.com/benjaminhgood/LTEE-metagenomic/issues/3)
replaced clang++ with g++ in make file.
in code/cpp folder run
$ make all
and moved resultant program "annotate_pvalues" to the base directory.

(10) MAPGD

Note from Will Shoemaker:
	You basically only need these three commands once you have your bam file sample.bam

	samtools sort -T /tmp/aln.sorted -o sample_sorted.bam sample.bam

	samtools view -H sample_sorted.bam > sample_sorted.header

	samtools mpileup -q 5 -Q 5 -B sample_sorted.bam \
	| mapgd proview -H sample_sorted.header | mapgd pool -a 20 -o sample_sorted.pol

	Though you can estimate allele frequencies for multiple samples at once to improve statistical power. I haven't done that but there are instructions on the mapgd README




#----------#
#  Credit  #
#----------#
The pipeline above is based on, and many scripts directly copied from:

Ben Good's repo
https://github.com/benjaminhgood/LTEE-metagenomic

Will Shoemaker's repos:
https://github.com/MURI2/Phylo_Evol_Timeseries
https://github.com/MURI2/Bacillus_Evol_Timeseries
