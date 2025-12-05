#!/usr/bin/env bash
#SBATCH --job-name=Nextflowpipe          # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-20:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=12               # Number of CPU cores per MPI task
#SBATCH --mem=256G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log


cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/E1_barcodeFollowUpTest/

#Load modules
module load minimap2/2.26-GCCcore-12.3.0
module load SAMtools/1.16.1-GCC-11.3.0

#Convert to .sam
minimap2 -ax map-ont /mnt/scratch/users/jps558/staging/hg38/hg38.fa JS_ES1_allpass_adapter_trimmed.fastq > aligned.sam
