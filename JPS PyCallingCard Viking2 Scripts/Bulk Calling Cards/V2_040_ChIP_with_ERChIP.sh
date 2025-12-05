#!/usr/bin/env bash
#SBATCH --job-name=ER_ChIP         # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-20:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=1               # Number of CPU cores per MPI task
#SBATCH --mem=50G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log

##Remove any single '#' at the end as these were put there to inactivate that line after successful completion of that line and later failure of a later line of code

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/EncodeChIP/

##module load and unload 
module purge
module load BEDOPS/2.4.41-foss-2021b
module load BEDTools/2.30.0-GCC-11.2.0

cd rep1/

##Convert to bed
bedtools bamtobed -i ENCFF365BIT.bam > ENCFF365BIT.bed
bedtools sort -i ENCFF365BIT.bed > ENCFF365BIT_sorted.bed

##Filter blacklist
bedops -n -1 ENCFF365BIT_sorted.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > ENCFF365BIT_filtered_sorted_hg38.bed

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/EncodeChIP/rep2/

##Convert to bed
bedtools bamtobed -i ENCFF063JMY.bam > ENCFF063JMY_hg38.bed
bedtools sort -i ENCFF063JMY_hg38.bed > ENCFF063JMY_sorted_hg38.bed

##Filter blacklist
bedops -n -1 ENCFF063JMY_sorted_hg38.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > ENCFF063JMY_filtered_sorted_hg38.bed


cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/EncodeChIP/input/

##Convert to bed
bedtools bamtobed -i ENCFF444JMD.bam > ENCFF444JMD_hg38.bed
bedtools sort -i ENCFF444JMD_hg38.bed > ENCFF444JMD_sorted_hg38.bed

##Filter blacklist
bedops -n -1 ENCFF444JMD_sorted_hg38.bed /mnt/scratch/users/jps558/staging/hg38/Blacklist/hg38-blacklist.v2.bed > ENCFF444JMD_filtered_sorted_hg38.bed


##Module load
module purge
module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/EncodeChIP/

#Call peaks for rep 1
macs2 callpeak -t rep1/ENCFF365BIT_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -q 0.05 \
    -n Encode_ER_ChIP_q0.05 \
    --outdir rep1/MACS2/BRD4_q0.05


macs2 callpeak -t rep1/ENCFF365BIT_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-9 \
    -n Encode_ER_rep1_ChIP_1e-9 \
    --outdir rep1/MACS2/1e-9

macs2 callpeak -t rep1/ENCFF365BIT_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-8 \
    -n Encode_ER_rep1_ChIP_1e-8 \
    --outdir rep1/MACS2/1e-8

macs2 callpeak -t rep1/ENCFF365BIT_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-7 \
    -n Encode_ER_rep1_ChIP_1e-7 \
    --outdir rep1/MACS2/1e-7

macs2 callpeak -t rep1/ENCFF365BIT_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-6 \
    -n Encode_ER_rep1_ChIP_1e-6 \
    --outdir rep1/MACS2/1e-6

macs2 callpeak -t rep1/ENCFF365BIT_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-5 \
    -n Encode_ER_rep1_ChIP_1e-5 \
    --outdir rep1/MACS2/1e-5

macs2 callpeak -t rep1/ENCFF365BIT_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-4 \
    -n Encode_ER_rep1_ChIP_1e-4 \
    --outdir rep1/MACS2/1e-4

#Call peaks for rep 2
macs2 callpeak -t rep2/ENCFF063JMY_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -q 0.05 \
    -n Encode_ER_ChIP_q0.05 \
    --outdir rep2/MACS2/BRD4_q0.05


macs2 callpeak -t rep2/ENCFF063JMY_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-9 \
    -n Encode_ER_rep2_ChIP_1e-9 \
    --outdir rep2/MACS2/1e-9

macs2 callpeak -t rep2/ENCFF063JMY_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-8 \
    -n Encode_ER_rep2_ChIP_1e-8 \
    --outdir rep2/MACS2/1e-8

macs2 callpeak -t rep2/ENCFF063JMY_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-7 \
    -n Encode_ER_rep2_ChIP_1e-7 \
    --outdir rep2/MACS2/1e-7

macs2 callpeak -t rep2/ENCFF063JMY_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-6 \
    -n Encode_ER_rep2_ChIP_1e-6 \
    --outdir rep2/MACS2/1e-6

macs2 callpeak -t rep2/ENCFF063JMY_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-5 \
    -n Encode_ER_rep2_ChIP_1e-5 \
    --outdir rep2/MACS2/1e-5

macs2 callpeak -t rep2/ENCFF063JMY_filtered_sorted_hg38.bed \
    -c input/ENCFF444JMD_filtered_sorted_hg38.bed \
    -f BED \
    -g hs \
    --keep-dup auto \
    -p 1e-4 \
    -n Encode_ER_rep2_ChIP_1e-4 \
    --outdir rep2/MACS2/1e-4









