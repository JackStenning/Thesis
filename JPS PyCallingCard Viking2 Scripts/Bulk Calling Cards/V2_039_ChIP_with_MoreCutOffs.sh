#!/usr/bin/env bash
#SBATCH --job-name=ChIP_cutoffs         # Job name
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

##Remove any single '#' at the end as these were put there to inactivate that line after successful completion of that line and later failure of a later line of code

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/MitraChIP

##Module load
module purge
module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

##Process SP1 FASTQ with Mitra control
macs2 callpeak -t ENCFF000PCT_filtered_sorted_hg38.bed \
    -c ENCFF000PBO_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-4 \
    -n SP1_ChIP_Mitra_input_1e-4 \
    --outdir MACS2/Mitra_input_1e-4

macs2 callpeak -t ENCFF000PCT_filtered_sorted_hg38.bed \
    -c ENCFF000PBO_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-5 \
    -n SP1_ChIP_Mitra_input_1e-5 \
    --outdir MACS2/Mitra_input_1e-5


macs2 callpeak -t ENCFF000PCT_filtered_sorted_hg38.bed \
    -c ENCFF000PBO_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-6 \
    -n SP1_ChIP_Mitra_input_1e-6 \
    --outdir MACS2/Mitra_input_1e-6

macs2 callpeak -t ENCFF000PCT_filtered_sorted_hg38.bed \
    -c ENCFF000PBO_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-7 \
    -n SP1_ChIP_Mitra_input_1e-7 \
    --outdir MACS2/Mitra_input_1e-7

macs2 callpeak -t ENCFF000PCT_filtered_sorted_hg38.bed \
    -c ENCFF000PBO_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-8 \
    -n SP1_ChIP_Mitra_input_1e-8 \
    --outdir MACS2/Mitra_input_1e-8

    
