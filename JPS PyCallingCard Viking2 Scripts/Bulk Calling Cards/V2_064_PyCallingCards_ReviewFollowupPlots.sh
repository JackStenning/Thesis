#!/bin/bash

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/

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

#Singal plot
from pycallingcards.plotting import signal_plot
import matplotlib.pyplot as plt

#Caluclate signal GEO

signalmtx = cc.pl.calculate_signal(peak_data_ER_WT2, 'SRR6652020_A90.bw')

#plot signal

signal_plot(
    signalmtx,
    fill_between=True,
    alpha=0.05,
    before=10000,
    after=10000,
    nbins=100,
    figsize=(8, 6),
    fontsize=10,
    color='red',
    textbelow=0,
    title='Log2(FC) Chip-seq Signal at High Stringency Peaks',
    save=False  # Ensure it doesn't auto-save, so we can modify it
)

# Get the current figure (since signal_plot() does not return fig)
fig = plt.gcf()  # Get the current figure
ax = plt.gca()   # Get the current axis

#Set Y-axis limits (Min: 1, Max: 4)
ax.set_ylim(1, 4)

# ? Customize Y-axis ticks
ax.tick_params(axis='y', labelsize=24)                 # Tick labels

#Add a solid red line at X = 0
ax.axvline(x=10000, color='red', linestyle='-', linewidth=1)  # Solid red vertical line

#Adjust X-axis label spacing (reduce gap)
ax.set_xlabel('Distance from Feature (bp)', labelpad=5)  # Reduce padding

#Save the modified figure
fig.savefig('Peakset2_GEO_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap

cc.pl.signal_heatmap(signalmtx, save = True)

###### ENCODE

os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_2')
#Caluclate signal Encode

signalmtx = cc.pl.calculate_signal(peak_data_ER_WT2, 'ENCFF063JMY.bw')

#plot signal

signal_plot(
    signalmtx,
    fill_between=True,
    alpha=0.05,
    before=10000,
    after=10000,
    nbins=100,
    figsize=(8, 6),
    fontsize=10,
    color='red',
    textbelow=0,
    title='Log2(FC) Chip-seq Signal at High Stringency Peaks',
    save=False  # Ensure it doesn't auto-save, so we can modify it
)

# Get the current figure (since signal_plot() does not return fig)
fig = plt.gcf()  # Get the current figure
ax = plt.gca()   # Get the current axis

#Set Y-axis limits (Min: 1, Max: 4)
ax.set_ylim(2, 8)

#Add a solid red line at X = 0
ax.axvline(x=10000, color='red', linestyle='-', linewidth=1)  # Solid red vertical line

#Adjust X-axis label spacing (reduce gap)
ax.set_xlabel('Distance from Feature (bp)', labelpad=5)  # Reduce padding

# ? Customize Y-axis ticks
ax.tick_params(axis='y', labelsize=24)                 # Tick labels

#Save the modified figure
fig.savefig('Peakset2_Encode_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap

cc.pl.signal_heatmap(signalmtx, colormap_vmax=8, save = True)





#Process Peakset 4
os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_4')

#Singal plot
#Caluclate signal

signalmtx4 = cc.pl.calculate_signal(peak_data_ER_WT4, 'SRR6652020_A90.bw')

#plot signal

signal_plot(
    signalmtx4,
    fill_between=True,
    alpha=0.05,
    before=10000,
    after=10000,
    nbins=100,
    figsize=(8, 6),
    fontsize=10,
    color='red',
    textbelow=0,
    title='Log2(FC) Chip-seq Signal at Low Stringency Peaks ',
    save=False  # Ensure it doesn't auto-save, so we can modify it
)

# Get the current figure (since signal_plot() does not return fig)
fig = plt.gcf()  # Get the current figure
ax = plt.gca()   # Get the current axis

#Set Y-axis limits (Min: 1, Max: 4)
ax.set_ylim(1, 4)

#Add a solid red line at X = 0
ax.axvline(x=10000, color='red', linestyle='-', linewidth=1)  # Solid red vertical line

#Adjust X-axis label spacing (reduce gap)
ax.set_xlabel('Distance from Feature (bp)', labelpad=5)  # Reduce padding

# ? Customize Y-axis ticks
ax.tick_params(axis='y', labelsize=24)                 # Tick labels

#Save the modified figure
fig.savefig('Peakset4_GEO_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)


#plot signal heatmap

cc.pl.signal_heatmap(signalmtx4, save = True)

###### ENCODE

os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_4')
#Caluclate signal Encode

signalmtx = cc.pl.calculate_signal(peak_data_ER_WT4, 'ENCFF063JMY.bw')

#plot signal

signal_plot(
    signalmtx,
    fill_between=True,
    alpha=0.05,
    before=10000,
    after=10000,
    nbins=100,
    figsize=(8, 6),
    fontsize=10,
    color='red',
    textbelow=0,
    title='Log2(FC) Chip-seq Signal at Low Stringency Peaks',
    save=False  # Ensure it doesn't auto-save, so we can modify it
)

# Get the current figure (since signal_plot() does not return fig)
fig = plt.gcf()  # Get the current figure
ax = plt.gca()   # Get the current axis

#Set Y-axis limits (Min: 1, Max: 4)
ax.set_ylim(2, 8)

#Add a solid red line at X = 0
ax.axvline(x=10000, color='red', linestyle='-', linewidth=1)  # Solid red vertical line

#Adjust X-axis label spacing (reduce gap)
ax.set_xlabel('Distance from Feature (bp)', labelpad=5)  # Reduce padding

# ? Customize Y-axis ticks
ax.tick_params(axis='y', labelsize=24)                 # Tick labels

#Save the modified figure
fig.savefig('Peakset4_Encode_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap

cc.pl.signal_heatmap(signalmtx, colormap_vmax=8, save = True)



# Draw Area
#Draw area test - 3 overlapping peaks - MAKE SURE THE PEAKS ARE RELEVENT -GREB 1 ETC
## MAKE SURE TO CHANGE FILE NAMES AS THEY ARE MADE TO ENSURE DATA IS NOT OVERWRITTEN ##

#Draw FOXA1 peak
cc.pl.draw_area("chr14", 37566280, 37567012, 50000, peak_data_ER_WT4, ER, "hg38", WT, color='green', font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 10000)
cc.pl.draw_area("chr14", 37566280, 37567012, 25000, peak_data_ER_WT4, ER, "hg38", WT, color='green', font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 5000)

#Draw GATA3 peak
cc.pl.draw_area("chr10", 8084674, 8085509, 10000, peak_data_ER_WT4, ER, "hg38", WT, color='purple', font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 5000)

#Draw GREB1 Peaks
cc.pl.draw_area("chr2", 11498705, 11499816, 50000, peak_data_ER_WT4, ER, "hg38", WT, color='red', font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 10000)
cc.pl.draw_area("chr2", 11654388, 11656833, 50000, peak_data_ER_WT4, ER, "hg38", WT, color='blue', font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 10000)



