#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=de-dup-pcr
cd  /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host 


/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_delta6-ANC-ANC-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-delta6-ANC-ANC-pl_S1_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-delta6-ANC-ANC-pl_S1_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_delta6-founder-T0-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-delta6-founder-T0-pl_S3_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-delta6-founder-T0-pl_S3_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_dSpoIIE-ANC-ANC-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-dSpoIIE-ANC-ANC-pl_S2_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-dSpoIIE-ANC-ANC-pl_S2_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_dSpoIIE-founder-T0-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-dSpoIIE-founder-T0-pl_S4_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-dSpoIIE-founder-T0-pl_S4_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNCt-L1-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L1-T14-pl_S20_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L1-T14-pl_S20_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNCt-L1-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L1-T14-rV_S50_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L1-T14-rV_S50_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNCt-L2-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L2-T14-pl_S21_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L2-T14-pl_S21_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNCt-L2-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L2-T14-rV_S51_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L2-T14-rV_S51_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNCt-L3-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L3-T14-pl_S22_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L3-T14-pl_S22_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNCt-L3-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L3-T14-rV_S52_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNCt-L3-T14-rV_S52_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNO-L1-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L1-T14-pl_S17_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L1-T14-pl_S17_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNO-L1-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L1-T14-rV_S47_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L1-T14-rV_S47_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNO-L2-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L2-T14-pl_S18_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L2-T14-pl_S18_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNO-L2-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L2-T14-rV_S48_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L2-T14-rV_S48_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNO-L3-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L3-T14-pl_S19_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L3-T14-pl_S19_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_SNO-L3-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L3-T14-rV_S49_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-SNO-L3-T14-rV_S49_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLCt-L1-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L1-T14-pl_S8_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L1-T14-pl_S8_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLCt-L1-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L1-T14-rS_S32_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L1-T14-rS_S32_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLCt-L1-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L1-T14-rV_S26_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L1-T14-rV_S26_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLCt-L2-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L2-T14-pl_S9_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L2-T14-pl_S9_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLCt-L2-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L2-T14-rS_S33_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L2-T14-rS_S33_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLCt-L2-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L2-T14-rV_S27_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L2-T14-rV_S27_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLCt-L3-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L3-T14-pl_S10_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L3-T14-pl_S10_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLCt-L3-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L3-T14-rS_S34_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L3-T14-rS_S34_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLCt-L3-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L3-T14-rV_S28_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLCt-L3-T14-rV_S28_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLO-L1-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L1-T14-pl_S5_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L1-T14-pl_S5_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLO-L1-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L1-T14-rS_S29_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L1-T14-rS_S29_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLO-L1-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L1-T14-rV_S23_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L1-T14-rV_S23_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLO-L2-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L2-T14-pl_S6_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L2-T14-pl_S6_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLO-L2-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L2-T14-rS_S30_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L2-T14-rS_S30_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLO-L2-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L2-T14-rV_S24_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L2-T14-rV_S24_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLO-L3-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L3-T14-pl_S7_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L3-T14-pl_S7_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLO-L3-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L3-T14-rS_S31_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L3-T14-rS_S31_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WLO-L3-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L3-T14-rV_S25_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WLO-L3-T14-rV_S25_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSCt-L1-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L1-T14-pl_S14_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L1-T14-pl_S14_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSCt-L1-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L1-T14-rS_S44_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L1-T14-rS_S44_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSCt-L1-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L1-T14-rV_S38_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L1-T14-rV_S38_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSCt-L2-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L2-T14-pl_S15_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L2-T14-pl_S15_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSCt-L2-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L2-T14-rS_S45_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L2-T14-rS_S45_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSCt-L2-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L2-T14-rV_S39_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L2-T14-rV_S39_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSCt-L3-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L3-T14-pl_S16_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L3-T14-pl_S16_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSCt-L3-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L3-T14-rS_S46_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L3-T14-rS_S46_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSCt-L3-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L3-T14-rV_S40_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSCt-L3-T14-rV_S40_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSO-L1-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L1-T14-pl_S11_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L1-T14-pl_S11_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSO-L1-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L1-T14-rS_S41_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L1-T14-rS_S41_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSO-L1-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L1-T14-rV_S35_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L1-T14-rV_S35_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSO-L2-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L2-T14-pl_S12_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L2-T14-pl_S12_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSO-L2-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L2-T14-rS_S42_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L2-T14-rS_S42_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSO-L2-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L2-T14-rV_S36_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L2-T14-rV_S36_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSO-L3-T14-pl.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L3-T14-pl_S13_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L3-T14-pl_S13_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSO-L3-T14-rS.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L3-T14-rS_S43_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L3-T14-rS_S43_R2_001.fastq -c 1 

/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq -i /N/slate/danschw/coevo-seedbank-seq/data/ddup-lists/host/input_list_WSO-L3-T14-rV.txt -t q -o /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L3-T14-rV_S37_R1_001.fastq -p /N/slate/danschw/coevo-seedbank-seq/data/ddup-fastq/host/ddup-GSF2865-WSO-L3-T14-rV_S37_R2_001.fastq -c 1 

