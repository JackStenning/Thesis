#!/usr/bin/env bash
#SBATCH --job-name=MitraCC_sort         # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-05:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=1               # Number of CPU cores per MPI task
#SBATCH --mem=50G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/MitraCC

#curl -L -o GSE148448_RAW.tar ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148448/suppl/GSE148448_RAW.tar
#tar -xvf GSE148448_RAW.tar

###Move data to the correct place
##Undirected PBase
#mkdir UndirectedPBase
#mv *HCT-116_PBase* UndirectedPBase
##SP1-PBase
#mkdir SP1-PBase
#mv *HCT-116_SP1-PBase* SP1-PBase

##module load and unload 
module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0
module load Miniconda3/23.5.2-0

#Unzip and sort
gunzip *.gz
for file in *ccf.txt; do
	bedtools sort -i "$file" > "$file"_sorted.qbed
	bedops -n -1 "$file"_sorted.qbed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > "$file"_filtered_sorted.qbed
	done
