#!/usr/bin/env bash
#SBATCH --job-name=ChIP_BRD         # Job name
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



##BRD4 Fastq, bw and Input fastq and bw
#FASTQ
#cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/MitraChIP/BRD4
#curl -L -o SRR2481799.fastq.gz https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR2481799&clipped=1

#BW - for both  SRR2481799 and SRR2481800 - downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73319
# Selected only the runs I wanted and manually transferred for ease

#input 
#mkdir input
#cd input
#curl -L -o SRR2481800.fastq.gz https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR2481800&clipped=1

##Load modules
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/1.20-GCC-13.2.0

#Process BRD4
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/MitraChIP/BRD4
gunzip *.gz
bwa mem -p -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa SRR2481799.fastq > SRR2481799_aligned.sam
samtools view -S -b SRR2481799_aligned.sam > SRR2481799.bam
samtools sort SRR2481799.bam -o SRR2481799_sorted_hg38.bam
samtools index SRR2481799_sorted_hg38.bam

#input
cd input
gunzip *.gz
bwa mem -p -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa SRR2481800.fastq > SRR2481800_aligned.sam
samtools view -S -b SRR2481800_aligned.sam > SRR2481800.bam
samtools sort SRR2481800.bam -o SRR2481800_sorted_hg38.bam
samtools index SRR2481800_sorted_hg38.bam

#Loading modules
module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0

#Process BRD4
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/MitraChIP/BRD4
bedtools bamtobed -i SRR2481799_sorted_hg38.bam > SRR2481799_hg38.bed
bedtools sort -i SRR2481799_hg38.bed > SRR2481799_sorted_hg38.bed
bedops -n -1 SRR2481799_sorted_hg38.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > SRR2481799_filtered_sorted_hg38.bed

#input
cd input
bedtools bamtobed -i SRR2481800_sorted_hg38.bam > SRR2481800_hg38.bed
bedtools sort -i SRR2481800_hg38.bed > SRR2481800_sorted_hg38.bed
bedops -n -1 SRR2481800_sorted_hg38.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > SRR2481800_filtered_sorted_hg38.bed

#Go back to /jack
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis