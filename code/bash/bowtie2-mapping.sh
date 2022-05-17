#!/bin/bash
#SBATCH --mail-user=williamrshoemaker@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=bowtie2-map
#SBATCH --mem=20G

# load dependencies
module load samtools
module load bowtie2/intel/2.4.2

module load gsl
module load htslib

#/geode/projects/iu/BL-BIO-Lennon-Lab/Data/0000_Schwartz/clean_coevo-seedbank-seq/coevo-seedbank-seq
# /N/slate/wrshoema
# index all the files
#wrshoema@carbonate.uits.iu.edu:/N/u/wrshoema/Carbonate/coevo-seedbank-seq/

# index fasta
#bowtie2-build /N/u/wrshoema/Carbonate/coevo-seedbank-seq/data/delta6-ANC.fa delta6-ANC
#bowtie2-build /N/u/wrshoema/Carbonate/coevo-seedbank-seq/data/SPO1-ANC.fa SPO1-ANC
#bowtie2-build -f /N/u/wrshoema/Carbonate/coevo-seedbank-seq/data/dspoIIE-ANC.fa dspoIIE-ANC


cd /N/u/wrshoema/Carbonate/coevo-seedbank-seq/data


# test
#sample_name=SNO-L3-T1
#fasta=SPO1-ANC
#R1=data/ddup-fastq/phage/ddup-GSF2842-9-SNO-L3-T1_S26_R1_001.fastq.gz
#R2=data/ddup-fastq/phage/ddup-GSF2842-9-SNO-L3-T1_S26_R2_001.fastq.gz

# create folder for bam output
while read sample; do
  #echo "$sample"
  unq_sample="$(echo "$sample" | cut -d ',' -f 1)"
  phageORhost="$(echo "$sample" | cut -d ',' -f 2)"
  seed_bank="$(echo "$sample" | cut -d ',' -f 7)"
  phage_trt="$(echo "$sample" | cut -d ',' -f 8)"
  subpop="$(echo "$sample" | cut -d ',' -f 9)"
  sample_name=${unq_sample}_${phageORhost}_${seed_bank}_${phage_trt}_${subpop}

  R1="$(echo "$sample" | cut -d ',' -f 13)"
  R2="$(echo "$sample" | cut -d ',' -f 14)"
  fasta="$(echo "$sample" | cut -d ',' -f 15)"
  fasta="$(echo "$fasta" | cut -d '/' -f 2 | cut -d '.' -f 1)"

  # skip header
  if [[ "$unq_sample" == 'unq.sample' ]]; then
    continue
  fi

  # revived_total, revived_spore, filtered_phage
  if [[ "$subpop" != 'revived_spore' ]]; then
    continue
  fi

  sam=/N/slate/wrshoema/coevo-seedbank-seq/data/sam/${sample_name}.sam
  bam=/N/slate/wrshoema/coevo-seedbank-seq/data/bam/${sample_name}.bam
  bam_sorted=/N/slate/wrshoema/coevo-seedbank-seq/data/bam/${sample_name}_sorted.bam
  bam_sorted_header=/N/slate/wrshoema/coevo-seedbank-seq/data/bam/${sample_name}_sorted.header
  pol=/N/slate/wrshoema/coevo-seedbank-seq/data/pol/${sample_name}

  #bowtie2 -x ${fasta} -1 /N/slate/wrshoema/coevo-seedbank-seq/${R1} -2 /N/slate/wrshoema/coevo-seedbank-seq/${R2} -S ${sam}

  #samtools view -S -b ${sam} > ${bam}
  #samtools sort -T /tmp/aln.sorted -o ${bam_sorted} ${bam}
  #samtools index ${bam_sorted}
  #samtools view -H ${bam_sorted} > ${bam_sorted_header}

  # try with -a 5 and see if it helps
  # old -a 15

  samtools mpileup -q 5 -Q 5 -B ${bam_sorted} | /N/u/wrshoema/Carbonate/MAPGD-master/bin/mapgd proview -H ${bam_sorted_header} | /N/u/wrshoema/Carbonate/MAPGD-master/bin/mapgd pool -a 5 -o ${pol}

  #rm ${sam}
  #rm ${bam}

done </N/u/wrshoema/Carbonate/coevo-seedbank-seq/data/samples_evolved.csv






# sbatch bowtie2-mapping.sh
