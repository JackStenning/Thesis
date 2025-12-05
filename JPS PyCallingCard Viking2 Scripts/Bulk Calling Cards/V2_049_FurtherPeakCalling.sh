#!/bin/bash

module load Apptainer/latest
module load Miniconda3/23.5.2-0

source activate CallingCardsPython

python3 -m pip install pyparsing twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy gnuplot
pip install "git+https://github.com/The-Mitra-Lab/pycallingcards.git" --upgrade
pip install git+https://github.com/cmatkhan/callingCardsTools.git

python
#From here code must be inputted into the terminal manually

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


peak_data_ER_WT3 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "peak_data_ER_WT3.bed")
peak_data_ER_WT3

peak_data_ER_WT4 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "peak_data_ER_WT4.bed")
peak_data_ER_WT4

peak_data_ER_WT5 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 100, record = True, save = "peak_data_ER_WT5.bed")
peak_data_ER_WT5

peak_data_ER_WT6 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 100, record = True, save = "peak_data_ER_WT6.bed")
peak_data_ER_WT6


peak_data_ER_WT7 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "peak_data_ER_WT7.bed")
peak_data_ER_WT7


peak_data_ER_WT8 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 1, record = True, save = "peak_data_ER_WT8.bed")
peak_data_ER_WT8

peak_data_ER_WT9 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 0.1, record = True, save = "peak_data_ER_WT9.bed")
peak_data_ER_WT9

peak_data_ER_WT10 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 0.1, record = True, save = "peak_data_ER_WT10.bed")
peak_data_ER_WT10

