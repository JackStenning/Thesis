#!/usr/bin/env bash
#SBATCH --job-name=Nextflowpipe          # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-10:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=12               # Number of CPU cores per MPI task
#SBATCH --mem=256G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/nextflow

srun --time=03:00:00 --partition=interactive --mem=10GB --cpus-per-task=4 --pty /bin/bash

module load Nextflow/23.04.2
module load Apptainer/latest


nextflow run epi2me-labs/wf-artic --help


#If you use epi2me-labs/wf-artic for your analysis please cite:

#* The nf-core framework
#  https://doi.org/10.1038/s41587-020-0439-x

#Update workflow
nextflow pull epi2me-labs/wf-single-cell

#Example code

#Below didnt work
nextflow run epi2me-labs/wf-single-cell \
  -r master \
  -profile singularity \
  --fastq 'JS_ES1_allpass_adapter_trimmed.fastq' \
  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
  --kit '3prime:v3' \
  --expected_cells 500 \
  --matrix_min_genes 100 \
  --matrix_min_cells 1 \
  --matrix_max_mito 30 \
  --out_dir 'result_full'

#Run on scCC trimmed samples
nextflow run epi2me-labs/wf-single-cell \
  -r master \
  -profile singularity \
  --fastq 'Extracted.fastq' \
  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
  --kit '3prime:v3' \
  --expected_cells 500 \
  --matrix_min_genes 100 \
  --matrix_min_cells 1 \
  --matrix_max_mito 30 \
  --out_dir 'result_pretrimmed'



nextflow run epi2me-labs/wf-single-cell \
  -profile singularity \
  --single_cell_sample_sheet my_primers_clean.tsv \
  --kit_config kit_config.csv \
  --fastq 'JS_ES1_allpass_adapter_trimmed.fastq' \
  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
  --out_dir 'result'