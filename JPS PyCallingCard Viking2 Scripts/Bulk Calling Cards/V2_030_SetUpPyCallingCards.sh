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
os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/ProcessPipe/ER')


ER_CGT = cc.rd.read_qbed("ER_test3CGT_1.qbed")
ER_CGT['group'] = 'ER_CGT'
ER_CGT

ER_CGT = cc.pp.clean_qbed(ER_CGT)

WT = cc.rd.read_qbed("WT_test3_1.qbed")
WT

WT = cc.pp.clean_qbed(WT)
WT

#Then try a number of different things - here are the notes for pvalue cut off args
#500-800 is good for step_size. pvalue_cutoffTTAA is the pvalue cutoff for TTAA data and pvalue_cutoffbg is pvalue cutoff for the background qbed data.
#pvalue_cutoffbg is recommended from 0.00001 to 0.01 and pvalue_cutoffTTAA is recommended from 0.001 to 0.1.
#The setting of pseudocounts is largely influenced by library size. For the first time of trial, it can be adjusted to 
((between 10^6 and 10^-5) x the number of insertions).


ERpeak_data_1 = cc.pp.call_peaks(ER_CGT, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1000, step_size = 500, pvalue_cutoffTTAA = 0.00001,
                             lam_win_size = 1000000, pseudocounts = 0.1, record = True, save = "ERpeak_1.bed")
ERpeak_data_1



#Above wrote nothing - so I am trying again with Juanru's code

ERpeak_data_2 = cc.pp.call_peaks(ER_CGT, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1000, step_size = 500, pvalue_cutoffTTAA = 0.00001,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "ERpeak_2.bed")
ERpeak_data_2

#Above also failed - I am confident its not an issue with the tool
