#!/usr/bin/env bash
#SBATCH --job-name=AndyChip         # Job name
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


cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/Andy_ChIP/

#Process ER
for file in *.fastq; do
	module purge
	module load BWA/0.7.17-GCCcore-11.3.0
	module load SAMtools/1.20-GCC-13.2.0
	bwa mem -p -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa "$file" > "$file"_aligned.sam
	samtools view -S -b "$file"_aligned.sam > "$file".bam
	samtools sort "$file".bam -o "$file"_sorted_hg38.bam
	samtools index "$file"_sorted_hg38.bam
	module purge
	module load BEDOPS/2.4.41-foss-2021b
	module load BEDTools/2.30.0-GCC-11.2.0
	bedtools bamtobed -i "$file"_sorted_hg38.bam > "$file".bed
	bedtools sort -i "$file".bed > "$file"_sorted.bed
	bedops -n -1 "$file"_sorted.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > "$file"_filtered_sorted.bed
	done

#Process Input
cd input

for file in *.fastq; do
	module purge
	module load BWA/0.7.17-GCCcore-11.3.0
	module load SAMtools/1.20-GCC-13.2.0
	bwa mem -p -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa "$file" > "$file"_aligned.sam
	samtools view -S -b "$file"_aligned.sam > "$file".bam
	samtools sort "$file".bam -o "$file"_sorted_hg38.bam
	samtools index "$file"_sorted_hg38.bam
	module purge
	module load BEDOPS/2.4.41-foss-2021b
	module load BEDTools/2.30.0-GCC-11.2.0
	bedtools bamtobed -i "$file"_sorted_hg38.bam > "$file".bed
	bedtools sort -i "$file".bed > "$file"_sorted.bed
	bedops -n -1 "$file"_sorted.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > "$file"_filtered_sorted.bed
	done
