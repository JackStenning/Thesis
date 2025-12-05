#!/usr/bin/env bash
#SBATCH --job-name=bamtobed         # Job name
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


cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/202306_JS_BulkCC/outputs/ER
module load BEDTools/2.31.0-GCC-12.3.0
module load BEDOPS/2.4.41-foss-2021b
bedtools bamtobed -i ER_HyPB_EKDL230008858_final.bam > ER_HyPB_EKDL230008858.bed
bedtools sort -i ER_HyPB_EKDL230008858.bed > ER_HyPB_EKDL230008858_sorted.bed
bedops -n -1 ER_HyPB_EKDL230008858_sorted.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > ER_HyPB_EKDL230008858_filtered_sorted.bed
cp ER_HyPB_EKDL230008858_filtered_sorted.bed /mnt/scratch/users/jps558/Stenning_Data_Analysis/


