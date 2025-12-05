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

#load modules
module load SRA-Toolkit/3.2.0-gompi-2024a

mkdir /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4_Zheng
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4_Zheng

#Get BRD4 input and ChIP fastq
# Download the .sra files
prefetch SRR16587391 SRR16587392

#91 is BRD, 92 is control

# Then convert to FASTQ
fasterq-dump SRR16587392 --split-files --threads 8
fasterq-dump SRR16587391 --split-files --threads 8

#Run for BRD4 ChIP
module purge
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/1.20-GCC-13.2.0
bwa mem -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa SRR16587391_1.fastq > SRR16587391_aligned.sam
samtools view -S -b SRR16587391_aligned.sam > SRR16587391.bam
samtools sort SRR16587391.bam -o SRR16587391_sorted_hg38.bam
samtools index SRR16587391_sorted_hg38.bam
module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0
bedtools bamtobed -i SRR16587391_sorted_hg38.bam > SRR16587391.bed
bedtools sort -i SRR16587391.bed > SRR16587391_sorted.bed
bedops -n -1 SRR16587391_sorted.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > SRR16587391_filtered_sorted.bed

#Run for control
module purge
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/1.20-GCC-13.2.0
bwa mem -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa SRR16587392_1.fastq > SRR16587392_aligned.sam
samtools view -S -b SRR16587392_aligned.sam > SRR16587392.bam
samtools sort SRR16587392.bam -o SRR16587392_sorted_hg38.bam
samtools index SRR16587392_sorted_hg38.bam
module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0
bedtools bamtobed -i SRR16587392_sorted_hg38.bam > SRR16587392.bed
bedtools sort -i SRR16587392.bed > SRR16587392_sorted.bed
bedops -n -1 SRR16587392_sorted.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > SRR16587392_filtered_sorted.bed




module purge
module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

macs2 callpeak -t SRR16587391_filtered_sorted.bed \
		-c SRR16587392_filtered_sorted.bed \
		-f BED \
		-g hs \
		--keep-dup auto \
		-q 0.05 \
		-n BRD4_control_Z_q0.05 \
		--outdir MACS2/Control_Z_q0.05



