#!/usr/bin/env bash
#SBATCH --job-name=BRD4Chip         # Job name
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

mkdir /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/Followup

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/Followup

## GET CHIAPET SIGNAL

#cd ~
#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#chmod +x Miniconda3-latest-Linux-x86_64.sh

#./Miniconda3-latest-Linux-x86_64.sh


#download interact files
wget -r -np -nH --cut-dirs=4 ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeGisChiaPet/ -P chia_pet/ -A "*Mcf7Eraa*.bigWig"

## CALCULATE BRD4 SIGNAL

module load BEDTools/2.31.0-GCC-12.3.0
module load deepTools/3.5.5

# Paths to BAM files
liu_bam="/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4/SRR14213436_sorted_hg38.bam"
zheng_bam="/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4_Zheng/SRR16587391_sorted_hg38.bam"

# Output directories (can be same as input)
liu_bw="/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4/SRR14213436_sorted_hg38.bw"
zheng_bw="/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4_Zheng/SRR16587391_sorted_hg38.bw"

# Run bamCoverage (from deepTools)
# --normalizeUsing CPM : counts per million mapped reads
# --binSize 10 : resolution of the coverage track (adjust if needed)
# --ignoreDuplicates : ignore PCR duplicates
# --extendReads : extend reads to fragment length (auto-detect)
# --numberOfProcessors : use multiple cores if available

bamCoverage -b "$liu_bam" -o "$liu_bw" --normalizeUsing CPM --binSize 10 --ignoreDuplicates --extendReads --numberOfProcessors 4
bamCoverage -b "$zheng_bam" -o "$zheng_bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    --ignoreDuplicates \
    --extendReads 200 \
    --numberOfProcessors 4

