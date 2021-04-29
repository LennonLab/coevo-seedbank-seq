#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=de-dup-pcr
cd  /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage 


/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L1-T1.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-7-SNO-L1-T1_S46_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-7-SNO-L1-T1_S46_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L1-T10.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-34-SNO-L1-T10_S58_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-34-SNO-L1-T10_S58_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L1-T14.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-8-SNO-L1-T14_S43_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-8-SNO-L1-T14_S43_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L1-T4.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-16-SNO-L1-T4_S29_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-16-SNO-L1-T4_S29_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L1-T7.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-25-SNO-L1-T7_S57_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-25-SNO-L1-T7_S57_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L2-T1.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-8-SNO-L2-T1_S47_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-8-SNO-L2-T1_S47_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L2-T10.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-35-SNO-L2-T10_S59_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-35-SNO-L2-T10_S59_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L2-T14.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-9-SNO-L2-T14_S44_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-9-SNO-L2-T14_S44_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L2-T4.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-17-SNO-L2-T4_S30_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-17-SNO-L2-T4_S30_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L2-T7.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-26-SNO-L2-T7_S35_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-26-SNO-L2-T7_S35_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L3-T1.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-9-SNO-L3-T1_S26_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-9-SNO-L3-T1_S26_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L3-T10.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-36-SNO-L3-T10_S43_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-36-SNO-L3-T10_S43_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L3-T14.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-10-SNO-L3-T14_S45_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-10-SNO-L3-T14_S45_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L3-T4.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-18-SNO-L3-T4_S31_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-18-SNO-L3-T4_S31_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SNO-L3-T7.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-27-SNO-L3-T7_S36_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-27-SNO-L3-T7_S36_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_SPO1-ANC-T0.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2742-1-SPO1-ANC-T0_S41_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2742-1-SPO1-ANC-T0_S41_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L1-T1.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-1-WLO-L1-T1_S44_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-1-WLO-L1-T1_S44_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L1-T10.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-28-WLO-L1-T10_S37_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-28-WLO-L1-T10_S37_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L1-T14.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-2-WLO-L1-T14_S38_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-2-WLO-L1-T14_S38_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L1-T4.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-10-WLO-L1-T4_S48_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-10-WLO-L1-T4_S48_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L1-T7.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-19-WLO-L1-T7_S54_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-19-WLO-L1-T7_S54_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L2-T1.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-2-WLO-L2-T1_S51_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-2-WLO-L2-T1_S51_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L2-T10.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-29-WLO-L2-T10_S38_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-29-WLO-L2-T10_S38_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L2-T14.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-3-WLO-L2-T14_S39_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-3-WLO-L2-T14_S39_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L2-T4.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-11-WLO-L2-T4_S27_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-11-WLO-L2-T4_S27_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L2-T7.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-20-WLO-L2-T7_S32_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-20-WLO-L2-T7_S32_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L3-T1.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-3-WLO-L3-T1_S52_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-3-WLO-L3-T1_S52_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L3-T10.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-30-WLO-L3-T10_S39_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-30-WLO-L3-T10_S39_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L3-T14.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2742-4-WLO-L3-T14_S42_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2742-4-WLO-L3-T14_S42_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L3-T4.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-12-WLO-L3-T4_S49_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-12-WLO-L3-T4_S49_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WLO-L3-T7.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-21-WLO-L3-T7_S55_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-21-WLO-L3-T7_S55_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L1-T1.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-4-WSO-L1-T1_S24_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-4-WSO-L1-T1_S24_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L1-T10.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-31-WSO-L1-T10_S40_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-31-WSO-L1-T10_S40_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L1-T14.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-5-WSO-L1-T14_S40_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-5-WSO-L1-T14_S40_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L1-T4.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-13-WSO-L1-T4_S28_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-13-WSO-L1-T4_S28_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L1-T7.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-22-WSO-L1-T7_S33_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-22-WSO-L1-T7_S33_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L2-T1.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-5-WSO-L2-T1_S45_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-5-WSO-L2-T1_S45_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L2-T10.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-32-WSO-L2-T10_S41_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-32-WSO-L2-T10_S41_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L2-T14.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-6-WSO-L2-T14_S41_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-6-WSO-L2-T14_S41_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L2-T4.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-14-WSO-L2-T4_S50_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-14-WSO-L2-T4_S50_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L2-T7.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-23-WSO-L2-T7_S56_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-23-WSO-L2-T7_S56_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L3-T1.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-6-WSO-L3-T1_S25_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-6-WSO-L3-T1_S25_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L3-T10.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-33-WSO-L3-T10_S42_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-33-WSO-L3-T10_S42_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L3-T14.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-7-WSO-L3-T14_S42_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2795-7-WSO-L3-T14_S42_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L3-T4.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-15-WSO-L3-T4_S53_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-15-WSO-L3-T4_S53_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/phage/input_list_WSO-L3-T7.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-24-WSO-L3-T7_S34_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/phage/ddup-GSF2842-24-WSO-L3-T7_S34_R2_001.fastq -c 1 

