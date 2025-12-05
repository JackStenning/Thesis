#!/usr/bin/env bash
#SBATCH --job-name=CCpipe          # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-10:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=12               # Number of CPU cores per MPI task
#SBATCH --mem=150G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log

module load Miniconda3/23.5.2-0
module load Nextflow/23.04.2
module load AlphaFold/2.3.1-foss-2022a
source activate CallingCardsPython3_8

#Setup
cd CC_Pipeline/
#mkdir Pipelineoutput
#callingcardstools barcode_table_to_json -t mammal_bc_table.tsv.tsv -r mammal

#Get Gtf
#curl -o hg38.refGene.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
#gunzip hg38.refGene.gtf.gz

nextflow run nf-core/callingcards \
    -profile default_mammals,conda \
    --input ESR1_test.csv \
    --fasta /mnt/scratch/users/jps558/staging/hg38/hg38.fa \
    --gtf hg38.refGene.gtf \
    --outdir Pipelineoutput

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis