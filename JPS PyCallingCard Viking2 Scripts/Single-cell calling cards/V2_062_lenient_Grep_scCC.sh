#!/usr/bin/env bash
#SBATCH --job-name=AndyGrep         # Job name
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


cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
#Find imperfect primer
for dir in *barcode*/; do
	cd "$dir"
	echo "$dir"
	for reads in *.fastq; do
	echo "$reads"
	grep -E "AAAGATAGTCTGCGTAAAATTGACGC.{0,3}$" "$reads" | grep -E "CACGACGCTCTTCCGATCT" > "${reads}AndyImperfectMatch.txt"
	cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
	done
done

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
#Find imperfect primer
for dir in *barcode*/; do
	cd "$dir"
	echo "$dir"
	for reads in *.fastq; do
	echo "$reads"
	grep -E "TTAACCCTAGAAAGATAGTCTGCGTAAAATTGACGC.{0,3}$" "$reads" | grep -E "CACGACGCTCTTCCGATCT" > "${reads}TTAAImperfectMatch.txt"
	cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
	done
done



