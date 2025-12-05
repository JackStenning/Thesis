#!/usr/bin/env bash
#SBATCH --job-name=AndyChip         # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-20:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=1               # Number of CPU cores per MPI task
#SBATCH --mem=150G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log


cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/Andy_ChIP/

#Get ER CHIP
#A90

curl -L -o SRR6652020_A90.fastq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR6652020&clipped=1

curl -L -o SRR6652018_A0.fastq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR6652018&clipped=1
gunzip *.gz

cd input

curl -L -o SRR6652038_A90.fastq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR6652038&clipped=1

curl -L -o SRR6652036_A0.fastq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR6652036&clipped=1

gunzip *.gz

cd ../
