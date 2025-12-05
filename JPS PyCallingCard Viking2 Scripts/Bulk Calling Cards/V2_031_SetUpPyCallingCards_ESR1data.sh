#!/bin/bash

module load Apptainer/latest
module load Miniconda3/23.5.2-0

source activate CallingCardsPython

python3 -m pip install pyparsing twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy gnuplot
pip install "git+https://github.com/The-Mitra-Lab/pycallingcards.git" --upgrade
pip install git+https://github.com/cmatkhan/callingCardsTools.git

python
import pycallingcards as cc
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
plt.rcParams['figure.dpi'] = 150

import os
os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year')


ER = cc.rd.read_qbed("ER_HyPB_EKDL230008858_sorted.qbed_filtered.qbed")
ER

ER = cc.pp.clean_qbed(ER)

WT = cc.rd.read_qbed("WT_HyPB_EKDL230008858_sorted.qbed_filtered.qbed")
WT

WT = cc.pp.clean_qbed(WT)
WT

peak_data_ER_WT1 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1000, step_size = 500, pvalue_cutoffTTAA = 0.00001,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "peak_data_ER_WT1.bed")
peak_data_ER_WT1

peak_data_ER_WT2 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "peak_data_ER_WT2.bed")
peak_data_ER_WT2

# call peaks for peak_data_ER1_HyPB
cc.tl.call_motif('peak_data_ER_WT1.bed', reference ="hg38",save_homer = "Homer/peak_data_ER_WT1", 
                 homer_path = "/home/juanru/miniconda3/bin/", num_cores=12)

# call peaks for peak_data_ER1_HyPB
cc.tl.call_motif('peak_data_ER_WT2.bed', reference ="hg38",save_homer = "Homer/peak_data_ER_WT2", 
                 homer_path = "/home/juanru/miniconda3/bin/", num_cores=12)

qbed = {"ER":ER, "WT": WT}
bed = {"peak_data_ER_WT1":peak_data_ER_WT1, "peak_data_ER_WT2":peak_data_ER_WT2}
cc.pl.WashU_browser_url(qbed = qbed,bed = bed,genome = 'hg38')

