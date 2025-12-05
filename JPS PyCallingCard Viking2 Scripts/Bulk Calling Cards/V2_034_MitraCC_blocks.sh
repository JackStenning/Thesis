#!/usr/bin/env bash
#SBATCH --job-name=MitraCC         # Job name
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

for file in *filtered_sorted.qbed; do
	python scripts/SegmentCCF.py "$file" | sed -e '/^\s*$/d' > "$file".blocks
	done
