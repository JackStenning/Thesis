#!/bin/bash
#SBATCH --job-name=scCCPyCC          # Job name
#SBATCH --partition=nodes                  # Partition
#SBATCH --time=0-25:00:00                  # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                         # Number of MPI tasks
#SBATCH --cpus-per-task=12                 # Number of CPU cores per task
#SBATCH --mem=256G                         # Total memory
#SBATCH --account=biol-cards-2023          # Project account
#SBATCH --mail-type=END,FAIL               # Mail events
#SBATCH --mail-user=jps558@york.ac.uk      # Email
#SBATCH --output=%x-%j.log                 # Standard output log
#SBATCH --error=%x-%j.err                  # Standard error log



#Load modules
module load Python/3.10.8-GCCcore-12.2.0
module load bzip2/1.0.8-GCCcore-14.2.0
module load HTSlib/1.19.1-GCC-13.2.0


cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/bedtools/
curl -L -o bedtools https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools
chmod +x bedtools
./bedtools --version


cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/E1_barcodeFollowUpTest

#Load environment
python -m venv ~/myenv
source ~/myenv/bin/activate

#Upgrade pip and install pure-Python dependencies
pip install --upgrade pip
python3 -m pip install pyparsing twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy gnuplot
pip install "git+https://github.com/The-Mitra-Lab/pycallingcards.git" --upgrade
pip install git+https://github.com/cmatkhan/callingCardsTools.git

python
#From here, the code must be inputted into the terminal manually
import shutil
import pycallingcards as cc
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
plt.rcParams['figure.dpi'] = 150

import os
os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/E1_barcodeFollowUpTest')
os.environ['PATH'] = '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/bedtools2:' + os.environ['PATH']


## Full results
Full_ER= cc.rd.read_qbed("result_full/JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf")
Full_ER

Full_ER= cc.pp.clean_qbed(Full_ER)


Full_ER_Peaks = cc.pp.call_peaks(Full_ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Full_ER_Peaks.bed")

Full_ER_Peaks

#UMI peaks
UMI_ER= cc.rd.read_qbed("result_full/JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf")
UMI_ER

UMI_ER= cc.pp.clean_qbed(UMI_ER)


UMI_ER_Peaks = cc.pp.call_peaks(UMI_ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./UMI_ER_Peaks.bed")

UMI_ER_Peaks

#Trim Peaks
Trim_ER= cc.rd.read_qbed("result_trim_minimum/JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf")
Trim_ER

Trim_ER= cc.pp.clean_qbed(Trim_ER)


Trim_ER_Peaks = cc.pp.call_peaks(Trim_ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Trim_ER_Peaks.bed")

Trim_ER_Peaks

#UMI peaks
Trim_UMI_ER= cc.rd.read_qbed("result_trim_minimum/JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf")
Trim_UMI_ER

Trim_UMI_ER= cc.pp.clean_qbed(Trim_UMI_ER)


Trim_UMI_ER_Peaks = cc.pp.call_peaks(Trim_UMI_ER, method = "MACCs", reference = "hg38", pvalue_cutoffbg = 0.01,
                             window_size = 1200, step_size = 800, pvalue_cutoffTTAA = 0.0005,
                             lam_win_size = 1000000, pseudocounts = 10, record = True, save = "./Trim_UMI_ER_Peaks.bed")

Trim_UMI_ER_Peaks

#### Full data analysis

# call Motif for Full data
cc.tl.call_motif('Full_ER_Peaks.bed', reference ="hg38",save_homer = "Homer/Full_ER_Peaks", 
                 homer_path = "/mnt/scratch/users/jps558/Stenning_Data_Analysis/homer/bin", num_cores=12)

# Combine Peak annotation
peak_annotationFull = cc.pp.annotation(Full_ER_Peaks, reference = "hg38", save_annotation = '/mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/E1_barcodeFollowUpTest/FullPeak_anno', bedtools_path= '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/bedtools/')
peak_annotationFull = cc.pp.combine_annotation(Full_ER_Peaks, peak_annotationFull)
peak_annotationFull
peak_annotationFull.to_csv('FullPeaks_annotated.csv', index=False)

#Draw GREB1 Peaks
cc.pl.draw_area("chr2", 11498705, 11499816, 50000, Full_ER_Peaks, Full_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'FullPeaks_chr2_11498705_11499816', bins = 500, example_length = 10000)
cc.pl.draw_area("chr2", 11654388, 11656833, 50000, Full_ER_Peaks, Full_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'FullPeaks_chr2_11654388_11656833', bins = 500, example_length = 10000)

#TFF1
cc.pl.draw_area("chr21", 42366284, 42366745, 50000, Full_ER_Peaks, Full_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'FullPeaks_chr21_42366284_42366745', bins = 500, example_length = 10000)

#KRT8/18
cc.pl.draw_area("chr12", 52897709, 52905392, 50000, Full_ER_Peaks, Full_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'FullPeaks_chr12_52897709_52905392', bins = 500, example_length = 10000, color = 'green')

#DNAJC3
cc.pl.draw_area("chr13", 95793873, 95794640, 50000, Full_ER_Peaks, Full_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'FullPeaks_chr13_95793873_95794640', bins = 500, example_length = 10000, color = 'blue')

#DDX5
cc.pl.draw_area("chr17", 64499547, 64506530, 50000, Full_ER_Peaks, Full_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'FullPeaks_chr17_64499547_64506530', bins = 500, example_length = 10000, color = 'blue')

#VMP1
cc.pl.draw_area("chr17", 59839695, 59840104, 50000, Full_ER_Peaks, Full_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'FullPeaks_chr17_59839695_59840104', bins = 500, example_length = 10000, color = 'blue')

#KCNJC3
cc.pl.draw_area("chr2", 154678774, 154679178, 50000, Full_ER_Peaks, Full_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'FullPeaks_chr2_154678774_154679178', bins = 500, example_length = 10000, color = 'green')

#Singal plot
from pycallingcards.plotting import signal_plot
import matplotlib.pyplot as plt

#Caluclate signal GEO

signalmtx = cc.pl.calculate_signal(Full_ER, '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/SRR6652020_A90.bw') 

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
    title='Log2(FC) Chip-seq Signal at Full ER peaks',
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
fig.savefig('FullResult_Peaks_GEO_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap

cc.pl.signal_heatmap(signalmtx, save = "FullResult_Peaks_GEO_heatmap.png")

#Caluclate signal Encode

signalmtx = cc.pl.calculate_signal(Full_ER, '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/ENCFF063JMY.bw') 

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
    title='Log2(FC) Chip-seq Signal at Full ER peaks',
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
fig.savefig('FullResult_Peaks_Encode_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap

cc.pl.signal_heatmap(signalmtx, save = "FullResult_Peaks_Encode_heatmap.png")


#### Full data analysis with UMI filter

# call Motif for UMI data
cc.tl.call_motif('UMI_ER_Peaks.bed', reference ="hg38",save_homer = "Homer/UMI_ER_Peaks", 
                 homer_path = "/mnt/scratch/users/jps558/Stenning_Data_Analysis/homer/bin", num_cores=12)

# Combine Peak annotation
peak_annotationUMI = cc.pp.annotation(UMI_ER_Peaks, reference = "hg38", save_annotation = '/mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/E1_barcodeFollowUpTest/UMIPeak_anno', bedtools_path= '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/bedtools2/')
peak_annotationUMI = cc.pp.combine_annotation(UMI_ER_Peaks, peak_annotationUMI)
peak_annotationUMI
peak_annotationUMI.to_csv('UMIPeaks_annotated.csv', index=False)

#Draw GREB1 Peaks
cc.pl.draw_area("chr2", 11498705, 11499816, 50000, UMI_ER_Peaks, UMI_ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 10000)
cc.pl.draw_area("chr2", 11654388, 11656833, 50000, UMI_ER_Peaks, UMI_ER, "hg38", WT, font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = True, bins = 500, example_length = 10000)

#TFF1
cc.pl.draw_area("chr21", 42366284, 42366745, 50000, UMI_ER_Peaks, UMI_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'UMIPeaks_chr21_42366284_42366745', bins = 500, example_length = 10000)

#KRT8/18
cc.pl.draw_area("chr12", 52897709, 52905392, 50000, UMI_ER_Peaks, UMI_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'UMIPeaks_chr12_52897709_52905392', bins = 500, example_length = 10000, color = 'green')

#DNAJC3
cc.pl.draw_area("chr13", 95793873, 95794640, 50000, UMI_ER_Peaks, UMI_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'UMIPeaks_chr13_95793873_95794640', bins = 500, example_length = 10000, color = 'blue')

#DDX5
cc.pl.draw_area("chr17", 64499547, 64506530, 50000, UMI_ER_Peaks, UMI_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'UMIPeaks_chr17_64499547_64506530', bins = 500, example_length = 10000, color = 'blue')

#VMP1
cc.pl.draw_area("chr17", 59839695, 59840104, 50000, UMI_ER_Peaks, UMI_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'UMIPeaks_chr17_59839695_59840104', bins = 500, example_length = 10000, color = 'blue')

#KCNJC3
cc.pl.draw_area("chr2", 154678774, 154679178, 50000, UMI_ER_Peaks, UMI_ER, "hg38", font_size=2, plotsize = [1,1,6], figsize = (30,12), peak_line = 2, save = 'UMIPeaks_chr2_154678774_154679178', bins = 500, example_length = 10000, color = 'green')


#Singal plot UMI 
from pycallingcards.plotting import signal_plot
import matplotlib.pyplot as plt

#Caluclate signal GEO

signalmtx = cc.pl.calculate_signal(UMI_ER, '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/SRR6652020_A90.bw') 

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
    title='Log2(FC) Chip-seq Signal at UMI ER peaks',
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
fig.savefig('UMIResult_Peaks_GEO_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap

cc.pl.signal_heatmap(signalmtx, save = "FullUMIResult_Peaks_GEO_heatmap.png")


#Caluclate signal Encode

signalmtx = cc.pl.calculate_signal(UMI_ER, '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/ENCFF063JMY.bw') 

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
    title='Log2(FC) Chip-seq Signal at UMI ER peaks',
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
fig.savefig('UMIResult_Peaks_Encode_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap

#plot signal heatmap

cc.pl.signal_heatmap(signalmtx, save = "FullUMIResult_Peaks_Encode_heatmap.png")





#Singal plot Trim
from pycallingcards.plotting import signal_plot
import matplotlib.pyplot as plt

#Caluclate signal GEO

signalmtx = cc.pl.calculate_signal(Trim_ER, '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/SRR6652020_A90.bw') 

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
    title='Log2(FC) Chip-seq Signal at Trim ER peaks',
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
fig.savefig('TrimResult_Peaks_GEO_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap

cc.pl.signal_heatmap(signalmtx, save = "TrimResult_Peaks_GEO_heatmap")


#Caluclate signal Encode

signalmtx = cc.pl.calculate_signal(Trim_ER, '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/ENCFF063JMY.bw') 

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
    title='Log2(FC) Chip-seq Signal at Trim ER peaks',
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
fig.savefig('TrimResult_Peaks_Encode_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap
cc.pl.signal_heatmap(signalmtx, save = "TrimResult_Peaks_Encode_heatmap")

#Singal plot Trim UMI
from pycallingcards.plotting import signal_plot
import matplotlib.pyplot as plt

#Caluclate signal GEO

signalmtx = cc.pl.calculate_signal(Trim_UMI_ER, '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/SRR6652020_A90.bw') 

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
    title='Log2(FC) Chip-seq Signal at Trim_UMI ER peaks',
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
fig.savefig('Trim_UMIResult_Peaks_GEO_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap

cc.pl.signal_heatmap(signalmtx, save = "TrimUMIResult_Peaks_GEO_heatmap")

#Caluclate signal Encode

signalmtx = cc.pl.calculate_signal(Trim_UMI_ER, '/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/ENCFF063JMY.bw') 

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
    title='Log2(FC) Chip-seq Signal at Trim_UMI ER peaks',
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
fig.savefig('Trim_UMIResult_Peaks_Encode_custom_signal_plot.png', dpi=300, bbox_inches='tight')

#Close the figure to free memory (important for batch processing)
plt.close(fig)

#plot signal heatmap
cc.pl.signal_heatmap(signalmtx, save = "TrimUMIResult_Peaks_Encode_heatmap")

