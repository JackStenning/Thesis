#!/bin/bash

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
os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year')
os.environ['PATH'] = '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/bedtools2:' + os.environ['PATH']


#Read Data into df

ER = cc.rd.read_qbed("ER_HyPB_EKDL230008858_sorted.qbed_filtered.qbed")
ER

ER = cc.pp.clean_qbed(ER)

WT = cc.rd.read_qbed("WT_HyPB_EKDL230008858_sorted.qbed_filtered.qbed")
WT

WT = cc.pp.clean_qbed(WT)
WT

#Call Peaks using optimal parameters

peak_data_ER_WT2 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.001,
                             window_size = 1200, step_size = 500, pvalue_cutoffTTAA = 0.00005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./peakset_2/peak_data_ER_WT2.bed")
peak_data_ER_WT2


peak_data_ER_WT4 = cc.pp.call_peaks(ER, WT, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./peakset_4/peak_data_ER_WT4.bed")

peak_data_ER_WT4


#Process Peakset 2
os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_2')

# call motif for peak_data_ER1_HyPB
cc.tl.call_motif('peak_data_ER_WT2.bed', reference ="hg38",save_homer = "Homer/peak_data_ER_WT2", 
                 homer_path = "/mnt/scratch/users/jps558/Stenning_Data_Analysis/homer/bin", num_cores=12)

# Combine Peak annotation
peak_annotation2 = cc.pp.annotation(peak_data_ER_WT2, reference = "hg38", save_annotation = '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_2/ER_WT2_anno', bedtools_path= '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/bedtools/')
peak_annotation2 = cc.pp.combine_annotation(peak_data_ER_WT2, peak_annotation2)
peak_annotation2
peak_annotation2.to_csv('Peakset2_annotated.csv', index=False)

#Draw area test - 3 overlapping peaks - MAKE SURE THE PEAKS ARE RELEVENT -GREB 1 ETC
## MAKE SURE TO CHANGE FILE NAMES AS THEY ARE MADE TO ENSURE DATA IS NOT OVERWRITTEN ##

cc.pl.draw_area("chr17", 50950089, 50951335, 100000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 50000)
cc.pl.draw_area("chr17", 50950089, 50951335, 20000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 5000)

#Draw FOXA1 peak
cc.pl.draw_area("chr14", 37566280, 37567012, 100000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 50000)
cc.pl.draw_area("chr14", 37566280, 37567012, 50000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 10000)

#Draw GATA3 peak
cc.pl.draw_area("chr10", 8084674, 8085509, 10000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 5000)

#Draw GREB1 Peaks
cc.pl.draw_area("chr2", 11498705, 11499816, 50000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 10000)
	
			
#Try whole peaks

cc.pl.whole_peaks(peak_data_ER_WT2, linewidth = 1, reference = "hg38", save = True)


#Singal plot
#Caluclate signal

mtx = cc.pl.calculate_signal(peak_data_ER_WT2, 'SRR6652020_A90.bw')

#plot signal

cc.pl.signal_plot(mtx, alpha = 0.2, save = True)

#plot signal heatmap

cc.pl.signal_heatmap(mtx, save = True)


#Process Peakset 4
os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_4')

# call peaks for peak_data_ER1_HyPB
cc.tl.call_motif('peak_data_ER_WT4.bed', reference ="hg38",save_homer = "Homer/peak_data_ER_WT4", 
                 homer_path = "/mnt/scratch/users/jps558/Stenning_Data_Analysis/homer/bin", num_cores=12)

# Combine Peak annotation
peak_annotation4 = cc.pp.annotation(peak_data_ER_WT4, reference = "hg38", save_annotation = '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_4/ER_WT4_anno', bedtools_path= '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/bedtools/')
peak_annotation4 = cc.pp.combine_annotation(peak_data_ER_WT4, peak_annotation4)
peak_annotation4
peak_annotation4.to_csv('Peakset4_annotated.csv', index=False)

#Draw area test - 3 overlapping peaks - MAKE SURE THE PEAKS ARE RELEVENT -GREB 1 ETC
## MAKE SURE TO CHANGE FILE NAMES AS THEY ARE MADE TO ENSURE DATA IS NOT OVERWRITTEN ##

cc.pl.draw_area("chr17", 50950089, 50951335, 100000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 50000)
cc.pl.draw_area("chr17", 50950089, 50951335, 20000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 5000)

#Draw FOXA1 peak
cc.pl.draw_area("chr14", 37566280, 37567012, 100000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 50000)
cc.pl.draw_area("chr14", 37566280, 37567012, 50000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 10000)

#Draw GATA3 peak
cc.pl.draw_area("chr10", 8084674, 8085509, 10000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 5000)

#Draw GREB1 Peaks
cc.pl.draw_area("chr2", 11498705, 11499816, 50000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 10000)
cc.pl.draw_area("chr2", 11654388, 11656833, 50000, peak_data_ER_WT4, ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 10000)

			
#Try whole peaks

cc.pl.whole_peaks(peak_data_ER_WT4, linewidth = 1, reference = "hg38", save = True)


#Singal plot
#Caluclate signal

mtx = cc.pl.calculate_signal(peak_data_ER_WT4, 'SRR6652020_A90.bw')

#plot signal

cc.pl.signal_plot(mtx, alpha = 0.2, save = True)

#plot signal heatmap

cc.pl.signal_heatmap(mtx, save = True)


