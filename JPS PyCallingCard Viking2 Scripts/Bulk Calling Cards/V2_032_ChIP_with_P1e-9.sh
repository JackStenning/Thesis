#!/usr/bin/env bash
#SBATCH --job-name=ChIP_1e-9         # Job name
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

##SP1 ChIP-Seq
curl  -L -o ENCFF000PCT.fastq.gz https://www.encodeproject.org/files/ENCFF000PCT/@@download/ENCFF000PCT.fastq.gz

##SP1 ChIP-Seq Mitra input
curl -L -o ENCFF000PBO.fastq.gz https://www.encodeproject.org/files/ENCFF000PBO/@@download/ENCFF000PBO.fastq.gz

##Load modules
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/1.20-GCC-13.2.0

##Unzip and align PCT
gunzip *.gz
bwa mem -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa ENCFF000PCT.fastq > ENCFF000PCT_aligned.sam

##Unzip and align PBO
gunzip *.gz
bwa mem -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa ENCFF000PBO.fastq > ENCFF000PBO_aligned.sam
  
##Convert to bam and sort
samtools view -S -b ENCFF000PCT_aligned.sam > ENCFF000PCT.bam
samtools sort ENCFF000PCT.bam -o ENCFF000PCT_sorted_hg38.bam
samtools index ENCFF000PCT_sorted_hg38.bam

  
##Convert to bam and sort
samtools view -S -b ENCFF000PBO_aligned.sam > ENCFF000PBO.bam
samtools sort ENCFF000PBO.bam -o ENCFF000PBO_sorted_hg38.bam
samtools index ENCFF000PBO_sorted_hg38.bam


##module load and unload 
module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0

##Convert to bed PCT
bedtools bamtobed -i ENCFF000PCT_sorted_hg38.bam > ENCFF000PCT_hg38.bed
bedtools sort -i ENCFF000PCT_hg38.bed > ENCFF000PCT_sorted_hg38.bed

##Filter blacklist
bedops -n -1 ENCFF000PCT_sorted_hg38.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > ENCFF000PCT_filtered_sorted_hg38.bed

##module load and unload 
module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0

##Convert to bed PBO
bedtools bamtobed -i ENCFF000PBO_sorted_hg38.bam > ENCFF000PBO_hg38.bed
bedtools sort -i ENCFF000PBO_hg38.bed > ENCFF000PBO_sorted_hg38.bed

##Filter blacklist
bedops -n -1 ENCFF000PBO_sorted_hg38.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > ENCFF000PBO_filtered_sorted_hg38.bed


##Module load
module purge
module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

##Process SP1 FASTQ with Mitra control
macs2 callpeak -t ENCFF000PCT_filtered_sorted_hg38.bed \
    -c ENCFF000PBO_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-9 \
    -n SP1_ChIP_Mitra_input_1e-9 \
    --outdir MACS2/Mitra_input_1e-9

macs2 callpeak -t ENCFF000PCT_filtered_sorted_hg38.bed \
    -c ENCFF000PBO_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -q 0.05 \
    -n SP1_ChIP_Mitra_input \
    --outdir MACS2/Mitra_input_asbefore

    
