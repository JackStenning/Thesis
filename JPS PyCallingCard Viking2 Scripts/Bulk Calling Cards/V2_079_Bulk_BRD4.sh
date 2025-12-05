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

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4

#curl -L -o GSE171908_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE171908&format=file"

#tar -xvf GSE171908_RAW.tar

#gunzip *.gz

#Get BRD4 input and ChIP fastq
# Download the .sra files
#prefetch SRR14213439 SRR14213436

# Then convert to FASTQ
#fasterq-dump SRR14213439 --split-files --threads 8
#fasterq-dump SRR14213436 --split-files --threads 8

#Run for BRD4 ChIP
module purge
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/1.20-GCC-13.2.0
bwa mem -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa SRR14213436_1.fastq SRR14213436_2.fastq > SRR14213436_aligned.sam
samtools view -S -b SRR14213436_aligned.sam > SRR14213436.bam
samtools sort SRR14213436.bam -o SRR14213436_sorted_hg38.bam
samtools index SRR14213436_sorted_hg38.bam
module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0
bedtools bamtobed -i SRR14213436_sorted_hg38.bam > SRR14213436.bed
bedtools sort -i SRR14213436.bed > SRR14213436_sorted.bed
bedops -n -1 SRR14213436_sorted.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > SRR14213436_filtered_sorted.bed

#Run for control
module purge
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/1.20-GCC-13.2.0
bwa mem -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa SRR14213439_1.fastq SRR14213439_2.fastq > SRR14213439_aligned.sam
samtools view -S -b SRR14213439_aligned.sam > SRR14213439.bam
samtools sort SRR14213439.bam -o SRR14213439_sorted_hg38.bam
samtools index SRR14213439_sorted_hg38.bam
module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0
bedtools bamtobed -i SRR14213439_sorted_hg38.bam > SRR14213439.bed
bedtools sort -i SRR14213439.bed > SRR14213439_sorted.bed
bedops -n -1 SRR14213439_sorted.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > SRR14213439_filtered_sorted.bed




module purge
module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

macs2 callpeak -t SRR14213436_filtered_sorted.bed \
		-c SRR14213439_filtered_sorted.bed \
		-f BED \
		-g hs \
		--keep-dup auto \
		-q 0.05 \
		-n BRD4_control_q0.05 \
		--outdir MACS2/Control_q0.05



