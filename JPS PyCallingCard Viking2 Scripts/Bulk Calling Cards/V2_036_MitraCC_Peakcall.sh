#!/usr/bin/env bash
#SBATCH --job-name=MitraCC_peaks         # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-20:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=12               # Number of CPU cores per MPI task
#SBATCH --mem=150G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/MitraCC

##module load and unload 
module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0
module load Apptainer/latest
module load Miniconda3/23.5.2-0

source activate CallingCardsPython
python3 -m pip install pyparsing twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy gnuplot

for dile in GSM4471638_HCT-116_HyPBase.ccf.txt_filtered_sorted.qbed; do
	python scripts/BBPeakCaller_BRD4.py -p 9 -d 12500 -i "$dile"_p9_intermediate.csv \
	"$dile" \
	"$dile".blocks \
	/mnt/scratch/users/jps558/staging/hg38/hg38_TTAA.bed \
	"$dile"_p9_peaks.bed
	done

#See how above runs before running this - find proper name for blocks
for mile in GSM4471639_HCT-116_SP1-HyPBase.ccf.txt_filtered_sorted.qbed; do
	python scripts/BBPeakCaller_TF.py -a 0.05 -m fdr_bh -d 250 -x 5000 \
	-i "$mile"_intermediate.csv \
	"$mile" \
	"$mile".blocks \
	GSM4471638_HCT-116_HyPBase.ccf.txt_filtered_sorted.qbed.blocks \
	"$mile"_fdr_bh_0.05_peaks.bed
	done

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis