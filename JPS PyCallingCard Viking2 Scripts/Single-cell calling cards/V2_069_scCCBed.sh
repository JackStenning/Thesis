#!/bin/bash

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/nextflow

#Load modules

module load BEDTools/2.31.0-GCC-12.3.0
module load BEDOPS/2.4.41-foss-2021b
module load SAMtools/1.16.1-GCC-11.3.0
module load Apptainer/latest
module load Python/3.11.3-GCCcore-12.3.0

## Full results

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/nextflow/result_full/JS_ES1_allpass_adapter_trimmed

python ../../AnnotateInsertionSites.py \
    --transposase PB \
   JS_ES1_allpass_adapter_trimmed.tagged.bam \
    /mnt/scratch/users/jps558/staging/hg38/hg38.2bit \
    Fullresult.CC_tagged.bam

cp JS_ES1_allpass_adapter_trimmed.tagged.bam.bai Fullresult.CC_tagged.bam.bai 

python ../../UMIFilter.py \
    -p 10x \
    -i Fullresult.CC_tagged.bam \
    --verbose \
    -o Fullresult.final.bam

cp Fullresult.CC_tagged.bam.bai Fullresult.final.bam.bai

python ../../BamToCallingCard.py \
    -b CB \
    -i Fullresult.CC_tagged.bam \
    > Fullresult.unsorted.ccf
grep -v '^chrM' Fullresult.unsorted.ccf > Fullresult_scCC_final.ccf
sort-bed Fullresult_scCC_final.ccf > Fullresult_scCC_final.sorted.ccf

python ../../BamToCallingCard.py \
    -b CB \
    -i Fullresult.final.bam \
    > Fullresult_UMIFilt.unsorted.ccf
sort-bed Fullresult_UMIFilt.unsorted.ccf > Fullresult_UMIFilt_scCC_final.ccf


## Pre-trimmed results with barcode

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/nextflow/result_pretrimmed_50Mit/Extracted

python ../../AnnotateInsertionSites.py \
    --transposase PB \
    Extracted.tagged.bam \
    /mnt/scratch/users/jps558/staging/hg38/hg38.2bit \
    Extracted.CC_tagged.bam

#Didnt work, only single UMI (maybe we can skip?)
#python ../../UMIFilter.py \
    -p 10x \
    -i Extracted.CC_tagged.bam \
    --verbose \
    -o Extracted.final.bam

cp Extracted.tagged.bam.bai Extracted.CC_tagged.bam.bai

python ../../BamToCallingCard.py \
    -b CB \
    -i Extracted.CC_tagged.bam \
    > Extracted.unsorted.ccf
sort-bed Extracted.unsorted.ccf > Extracted_scCC_final.ccf


## Pre-trimmed results with barcode (min 1)

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/nextflow/result_pretrimmed_min1/Extracted_min1

python ../../AnnotateInsertionSites.py \
    --transposase PB \
    Extracted_min1.tagged.bam \
    /mnt/scratch/users/jps558/staging/hg38/hg38.2bit \
    Extracted_min1.CC_tagged.bam

#Didnt work, only single UMI (maybe we can skip?)
#python ../../UMIFilter.py \
    -p 10x \
    -i Extracted_min1.CC_tagged.bam \
    --verbose \
    -o Extracted_min1.final.bam 

cp Extracted_min1.tagged.bam.bai Extracted_min1.CC_tagged.bam.bai

python ../../BamToCallingCard.py \
    -b CB \
    -i Extracted_min1.CC_tagged.bam \
    > Extracted_min1.unsorted.ccf
sort-bed Extracted_min1.unsorted.ccf > Extracted_min1_scCC_final.ccf


## Pre-trimmed results without barcode full result

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/nextflow/result_nobc/nobcExtracted_Final

python ../../AnnotateInsertionSites.py \
    --transposase PB \
    nobcExtracted_Final.tagged.bam \
    /mnt/scratch/users/jps558/staging/hg38/hg38.2bit \
    nobcExtracted.CC_tagged.bam

#Didnt work, only single UMI (maybe we can skip?)
#python ../../UMIFilter.py \
    -p 10x \
    -i nobcExtracted.CC_tagged.bam \
    --verbose \
    -o nobcExtracted.final.bam

cp nobcExtracted_Final.tagged.bam.bai nobcExtracted.CC_tagged.bam.bai

python ../../BamToCallingCard.py \
    -b CB \
    -i nobcExtracted.CC_tagged.bam \
    > nobcExtracted.unsorted.ccf
sort-bed nobcExtracted.unsorted.ccf > nobcExtracted_scCC_final.ccf

## Pre-trimmed results without barcode minimal epi2me

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/nextflow/result_nobc_min/nobcExtracted_Final

python ../../AnnotateInsertionSites.py \
    --transposase PB \
    nobcExtracted_Final.tagged.bam \
    /mnt/scratch/users/jps558/staging/hg38/hg38.2bit \
    nobcExtracted_min.CC_tagged.bam

#Didnt work, only single UMI (maybe we can skip?)
#python ../../UMIFilter.py \
    -p 10x \
    -i nobcExtracted_min.CC_tagged.bam \
    --verbose \
    -o nobcExtracted_min.final.bam

cp nobcExtracted_Final.tagged.bam.bai nobcExtracted_min.CC_tagged.bam.bai

python ../../BamToCallingCard.py \
    -b CB \
    -i nobcExtracted_min.CC_tagged.bam \
    > nobcExtracted_min.unsorted.ccf
sort-bed nobcExtracted_min.unsorted.ccf > nobcExtracted_min_scCC_final.ccf



#Set up Python
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
os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/nextflow/FinalCCF')

## Full results
ER = cc.rd.read_qbed("Fullresult_scCC_final.sorted.ccf")
ER

ER = cc.pp.clean_qbed(ER)

peak_data_ER_WT2 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Fullresult/peak_data_ER_test2.bed")
peak_data_ER_WT2

peak_data_ER_WT4 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Fullresult/peak_data_ER_test4.bed")

peak_data_ER_WT4

#With UMI filter

ER = cc.rd.read_qbed("Fullresult_UMIFilt_scCC_final.ccf")
ER

ER = cc.pp.clean_qbed(ER)

peak_data_ER_WT2 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./FullresultUMI/peak_data_ER_test2.bed")
peak_data_ER_WT2

peak_data_ER_WT4 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./FullresultUMI/peak_data_ER_test4.bed")

peak_data_ER_WT4



#### Pre-trimmed results with barcode
ER = cc.rd.read_qbed("Extracted_scCC_final.ccf")
ER

ER = cc.pp.clean_qbed(ER)

peak_data_ER_WT2 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Extracted_scCC/peak_data_ER_test2.bed")
peak_data_ER_WT2

peak_data_ER_WT4 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Extracted_scCC/peak_data_ER_test4.bed")

peak_data_ER_WT4

#### Pre-trimmed results with barcode minimum 1 cell

ER = cc.rd.read_qbed("Extracted_min1_scCC_final.ccf")
ER

ER = cc.pp.clean_qbed(ER)

peak_data_ER_WT2 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Extracted_min1/peak_data_ER_test2.bed")
peak_data_ER_WT2

peak_data_ER_WT4 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Extracted_min1/peak_data_ER_test4.bed")

peak_data_ER_WT4

#### Pre-trimmed results without barcode

ER = cc.rd.read_qbed("nobcExtracted_scCC_final.ccf")
ER

ER = cc.pp.clean_qbed(ER)

peak_data_ER_WT2 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./nobc/peak_data_ER_test2.bed")
peak_data_ER_WT2

peak_data_ER_WT4 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./nobc/peak_data_ER_test4.bed")

peak_data_ER_WT4

#### Pre-trimmed results without barcode min1

ER = cc.rd.read_qbed("nobcExtracted_min_scCC_final.ccf")
ER

ER = cc.pp.clean_qbed(ER)

peak_data_ER_WT2 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./nobcMin/peak_data_ER_test2.bed")
peak_data_ER_WT2

peak_data_ER_WT4 = cc.pp.call_peaks(ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./nobcMin/peak_data_ER_test4.bed")

peak_data_ER_WT4


