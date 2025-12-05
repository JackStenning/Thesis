#!/usr/bin/env bash
#SBATCH --job-name=ChIP_BRD_1e-9         # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-02:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=1               # Number of CPU cores per MPI task
#SBATCH --mem=50G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log

##Remove any single '#' at the end as these were put there to inactivate that line after successful completion of that line and later failure of a later line of code

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/MitraChIP/BRD4

##Module load
module purge
module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

##Process SP1 FASTQ with Mitra control
macs2 callpeak -t SRR2481799_filtered_sorted_hg38.bed \
    -c input/SRR2481800_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-9 \
    -n BRD4_ChIP_1e-9 \
    --outdir MACS2/BRD4_1e-9

macs2 callpeak -t SRR2481799_filtered_sorted_hg38.bed \
    -c input/SRR2481800_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -q 0.05 \
    -n BRD4_ChIP_q0.05 \
    --outdir MACS2/BRD4_q0.05

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis
    
