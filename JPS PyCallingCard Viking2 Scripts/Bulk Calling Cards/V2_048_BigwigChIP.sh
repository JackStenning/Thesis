#!/usr/bin/env bash
#SBATCH --job-name=BedtoBW         # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-20:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=3               # Number of CPU cores per MPI task
#SBATCH --mem=100G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log

module load BEDTools/2.31.0-GCC-12.3.0
module load BEDOPS/2.4.41-foss-2021b

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/Andy_ChIP/

bedtools genomecov -bg -i SRR6652020_A90.fastq_filtered_sorted.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38_genomeFile.txt > SRR6652020_A90.bedgraph
module load Kent_tools/442-GCC-11.3.0
sort -k1,1 -k2,2n SRR6652020_A90.bedgraph > SRR6652020_A90_sorted.bedgraph
bedGraphToBigWig SRR6652020_A90_sorted.bedgraph /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes SRR6652020_A90.bw

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/EncodeChIP/rep1

bedtools genomecov -bg -i ENCFF365BIT_filtered_sorted_hg38.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38_genomeFile.txt > ENCFF365BIT.bedgraph
module load Kent_tools/442-GCC-11.3.0
sort -k1,1 -k2,2n ENCFF365BIT.bedgraph > ENCFF365BIT_sorted.bedgraph
bedGraphToBigWig ENCFF365BIT_sorted.bedgraph /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes ENCFF365BIT.bw

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/EncodeChIP/rep2

bedtools genomecov -bg -i ENCFF063JMY_filtered_sorted_hg38.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38_genomeFile.txt > ENCFF063JMY.bedgraph
module load Kent_tools/442-GCC-11.3.0
sort -k1,1 -k2,2n ENCFF063JMY.bedgraph > ENCFF063JMY_sorted.bedgraph
bedGraphToBigWig ENCFF063JMY_sorted.bedgraph /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes ENCFF365BIT.bw
