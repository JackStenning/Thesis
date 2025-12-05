#!/usr/bin/env bash
#SBATCH --job-name=Window         # Job name
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

module load BEDTools/2.30.0-GCC-11.2.0

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/

#Define variables
root=/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/
mkdir results/
result=/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/results/

for dir in *ChIP/; do
	cd "$dir"
	for mac in MACS2*/; do
		cd "$mac"
		for pval in */; do
			cd "$pval"
			rm *p9_peaks.bed
			rm *p30_peaks.bed
			rm *p62_peaks.bed
			rm *peaks_r.bed
			cp "$root"*peak* .
			for file in *.narrowPeak; do
			  cp "$file" "$root""$dir""$mac""$pval""$(basename "$file" .narrowPeak)NP.bed"
			done
			rm *trim*
			rm narrowpeak.bed			
			for CC in peak_data*; do
				for chip in *NP.bed; do
					bedtools window -w 1000 -c -a "$CC" -b "$chip" > "$result"1kbp_A_"$CC"_B_"$chip"
					done
				done

			cd ../
			done
		cd ../
		done
	cd ../
	done
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/