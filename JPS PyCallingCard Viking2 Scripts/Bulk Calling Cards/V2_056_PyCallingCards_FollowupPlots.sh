#!/bin/bash

module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0
module load Apptainer/latest
module load Miniconda3/23.5.2-0

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/202306_JS_BulkCC/outputs/ER

#Combine files form flow cells
cat ER1_EKDL230008858-1A_HF552DSX7_L2_1_.qbed ER1_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combined/ER1_EKDL230008858.qbed
cat ER2_EKDL230008858-1A_HF552DSX7_L2_1_.qbed ER2_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combined/ER2_EKDL230008858.qbed
cat ER3_EKDL230008858-1A_HF552DSX7_L2_1_.qbed ER3_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combined/ER3_EKDL230008858.qbed
cat ER4_EKDL230008858-1A_HF552DSX7_L2_1_.qbed ER4_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combined/ER4_EKDL230008858.qbed
cat ER5_EKDL230008858-1A_HF552DSX7_L2_1_.qbed ER5_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combined/ER5_EKDL230008858.qbed
cat ER6_EKDL230008858-1A_HF552DSX7_L2_1_.qbed ER6_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combined/ER6_EKDL230008858.qbed

#Combine replicates together
cat combined/ER1_EKDL230008858.qbed combined/ER2_EKDL230008858.qbed combined/ER3_EKDL230008858.qbed > combined/ER1-3_EKDL230008858.qbed
cat combined/ER5_EKDL230008858.qbed combined/ER6_EKDL230008858.qbed > combined/ER5-6_EKDL230008858.qbed
cat combined/ER1_EKDL230008858.qbed combined/ER2_EKDL230008858.qbed combined/ER3_EKDL230008858.qbed combined/ER5_EKDL230008858.qbed combined/ER6_EKDL230008858.qbed > combined/ER_HyPB_EKDL230008858.qbed


cp -R combined/ /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/Replicate_analysis/

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/202306_JS_BulkCC/outputs/WT

mkdir combinedWT/ 

#Combine files form flow cells
cat WT1_EKDL230008858-1A_HF552DSX7_L2_1_.qbed WT1_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combinedWT/WT1_EKDL230008858.qbed
cat WT2_EKDL230008858-1A_HF552DSX7_L2_1_.qbed WT2_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combinedWT/WT2_EKDL230008858.qbed
cat WT3_EKDL230008858-1A_HF552DSX7_L2_1_.qbed WT3_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combinedWT/WT3_EKDL230008858.qbed
cat WT4_EKDL230008858-1A_HF552DSX7_L2_1_.qbed WT4_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combinedWT/WT4_EKDL230008858.qbed
cat WT5_EKDL230008858-1A_HF552DSX7_L2_1_.qbed WT5_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combinedWT/WT5_EKDL230008858.qbed
cat WT6_EKDL230008858-1A_HF552DSX7_L2_1_.qbed WT5_EKDL230008858-1A_HFK3FDSX7_L2_1_.qbed > combinedWT/WT6_EKDL230008858.qbed


#Combine files replicates together
cat combinedWT/WT1_EKDL230008858.qbed combinedWT/WT2_EKDL230008858.qbed combinedWT/WT3_EKDL230008858.qbed > combinedWT/WT1-3_EKDL230008858.qbed
cat combinedWT/WT5_EKDL230008858.qbed combinedWT/WT6_EKDL230008858.qbed > combinedWT/WT5-6_EKDL230008858.qbed
cat combinedWT/WT1_EKDL230008858.qbed combinedWT/WT2_EKDL230008858.qbed combinedWT/WT3_EKDL230008858.qbed combinedWT/WT5_EKDL230008858.qbed combinedWT/WT6_EKDL230008858.qbed> combinedWT/WT_HyPB_EKDL230008858.qbed



cp -R combinedWT/ /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/Replicate_analysis/


#Sort .qbed
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/Replicate_analysis/

for dir in combined*/; do
	cd $dir
	for bed in *.qbed; do
	bedtools sort -i $bed > sorted_$bed
	bedtools subtract -A -a  sorted_$bed -b /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > filtered_sorted$bed
	done
	cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/Replicate_analysis/
done


#peakcall
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis

module load Apptainer/latest
module load Miniconda3/23.5.2-0

source activate CallingCardsPython

python3 -m pip install pyparsing twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy gnuplot
pip install "git+https://github.com/The-Mitra-Lab/pycallingcards.git" --upgrade
pip install git+https://github.com/cmatkhan/callingCardsTools.git

python
#From here code must be inputted into the terminal manually
import shutil
import pycallingcards as cc
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
plt.rcParams['figure.dpi'] = 150

import os
os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/Replicate_analysis')

#Read combined Data into df

ERcombi = cc.rd.read_qbed("combined/filtered_sortedER_HyPB_EKDL230008858.qbed")
ERcombi

ERcombi = cc.pp.clean_qbed(ERcombi)

WTcombi = cc.rd.read_qbed("combinedWT/filtered_sortedWT_HyPB_EKDL230008858.qbed")
WTcombi

WTcombi = cc.pp.clean_qbed(WTcombi)
WTcombi



#Call Peaks using optimal parameters

peak_data_ER_WTcombi2 = cc.pp.call_peaks(ERcombi, WTcombi, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTcombi_2.bed")
peak_data_ER_WTcombi2


peak_data_ER_WTcombi4 = cc.pp.call_peaks(ERcombi, WTcombi, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTcombi_4.bed")

peak_data_ER_WTcombi4



#Read Rep1 Data into df

ER1 = cc.rd.read_qbed("combined/filtered_sortedER1_EKDL230008858.qbed")
ER1

ER1 = cc.pp.clean_qbed(ER1)

WT1 = cc.rd.read_qbed("combinedWT/filtered_sortedWT1_EKDL230008858.qbed")
WT1

WT1 = cc.pp.clean_qbed(WT1)
WT1



#Call Peaks using optimal parameters

peak_data_ER_WT1_2 = cc.pp.call_peaks(ER1, WT1, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep1_2.bed")
peak_data_ER_WT1_2


peak_data_ER_WT1_4 = cc.pp.call_peaks(ER1, WT1, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep1_4.bed")

peak_data_ER_WT1_4


#Read Rep2 data into df

ER2 = cc.rd.read_qbed("combined/filtered_sortedER2_EKDL230008858.qbed")
ER2

ER2 = cc.pp.clean_qbed(ER2)

WT2 = cc.rd.read_qbed("combinedWT/filtered_sortedWT2_EKDL230008858.qbed")
WT2

WT2 = cc.pp.clean_qbed(WT2)
WT2



#Call Peaks using optimal parameters

peak_data_ER_WT2_2 = cc.pp.call_peaks(ER2, WT2, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep2_2.bed")
peak_data_ER_WT2_2


peak_data_ER_WT2_4 = cc.pp.call_peaks(ER2, WT2, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep2_4.bed")

peak_data_ER_WT2_4


#Read Rep3 data into df

ER3 = cc.rd.read_qbed("combined/filtered_sortedER3_EKDL230008858.qbed")
ER3

ER3 = cc.pp.clean_qbed(ER3)

WT3 = cc.rd.read_qbed("combinedWT/filtered_sortedWT3_EKDL230008858.qbed")
WT3

WT3 = cc.pp.clean_qbed(WT3)
WT3



#Call Peaks using optimal parameters

peak_data_ER_WT3_2 = cc.pp.call_peaks(ER3, WT3, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep3_2.bed")
peak_data_ER_WT3_2


peak_data_ER_WT3_4 = cc.pp.call_peaks(ER3, WT3, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep3_4.bed")

peak_data_ER_WT3_4



#Read Rep5 data into df

ER5 = cc.rd.read_qbed("combined/filtered_sortedER5_EKDL230008858.qbed")
ER5

ER5 = cc.pp.clean_qbed(ER5)

WT5 = cc.rd.read_qbed("combinedWT/filtered_sortedWT5_EKDL230008858.qbed")
WT5

WT5 = cc.pp.clean_qbed(WT5)
WT5



#Call Peaks using optimal parameters

peak_data_ER_WT5_2 = cc.pp.call_peaks(ER5, WT5, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep5_2.bed")
peak_data_ER_WT5_2


peak_data_ER_WT5_4 = cc.pp.call_peaks(ER5, WT5, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep5_4.bed")

peak_data_ER_WT5_4



#Read Rep6 data into df

ER6 = cc.rd.read_qbed("combined/filtered_sortedER6_EKDL230008858.qbed")
ER6

ER6 = cc.pp.clean_qbed(ER6)

WT6 = cc.rd.read_qbed("combinedWT/filtered_sortedWT6_EKDL230008858.qbed")
WT6

WT6 = cc.pp.clean_qbed(WT6)
WT6



#Call Peaks using optimal parameters

peak_data_ER_WT6_2 = cc.pp.call_peaks(ER6, WT6, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep6_2.bed")
peak_data_ER_WT6_2


peak_data_ER_WT6_4 = cc.pp.call_peaks(ER6, WT6, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep6_4.bed")

peak_data_ER_WT6_4


#Read Rep1-3 data into df

ER1_3 = cc.rd.read_qbed("combined/filtered_sortedER1-3_EKDL230008858.qbed")
ER1_3

ER1_3 = cc.pp.clean_qbed(ER1_3)

WT1_3 = cc.rd.read_qbed("combinedWT/filtered_sortedWT1-3_EKDL230008858.qbed")
WT1_3

WT1_3 = cc.pp.clean_qbed(WT1_3)
WT1_3



#Call Peaks using optimal parameters

peak_data_ER_WT1_3_2 = cc.pp.call_peaks(ER1_3, WT1_3, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep1_3_2.bed")
peak_data_ER_WT1_3_2


peak_data_ER_WT1_3_4 = cc.pp.call_peaks(ER1_3, WT1_3, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep1_3_4.bed")

peak_data_ER_WT1_3_4


#Read Rep5-6 data into df

ER5_6 = cc.rd.read_qbed("combined/filtered_sortedER5-6_EKDL230008858.qbed")
ER5_6

ER5_6 = cc.pp.clean_qbed(ER5_6)

WT5_6 = cc.rd.read_qbed("combinedWT/filtered_sortedWT5-6_EKDL230008858.qbed")
WT5_6

WT5_6 = cc.pp.clean_qbed(WT5_6)
WT5_6



#Call Peaks using optimal parameters

peak_data_ER_WT5_6_2 = cc.pp.call_peaks(ER5_6, WT5_6, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep5_6_2.bed")
peak_data_ER_WT5_6_2


peak_data_ER_WT5_6_4 = cc.pp.call_peaks(ER5_6, WT5_6, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Peakcall/peak_data_ER_WTrep5_6_4.bed")

peak_data_ER_WT5_6_4


