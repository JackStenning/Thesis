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

############################# THIS SCRIPT IS THE SAME AS V2_068_1 BUT HAS CODE THAT WAS RUN ON A PRIOR DATE REMOVED. WHEN REPLICATING ANALYSIS, ONLY V2_068_1 IS NEEDED ########################

module load Nextflow/23.04.2
module load Apptainer/latest


nextflow run epi2me-labs/wf-artic --help


#If you use epi2me-labs/wf-artic for your analysis please cite:

#* The nf-core framework
#  https://doi.org/10.1038/s41587-020-0439-x

#Update workflow
nextflow pull epi2me-labs/wf-single-cell -r master

# Run on scCC trimmed samples with barcode (minus code already run in V2_068_1...)

nextflow run epi2me-labs/wf-single-cell \
  -r master \
  -profile singularity \
  --fastq 'Extracted_min1.fastq' \
  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
  --kit '3prime:v3' \
  --expected_cells 1 \
  --matrix_min_genes 1 \
  --matrix_min_cells 1 \
  --matrix_max_mito 50 \
  --out_dir 'result_pretrimmed_min1'


#Run on trimmmed without barcode
#Full
nextflow run epi2me-labs/wf-single-cell \
  -r master \
  -profile singularity \
  --fastq 'nobcExtracted_Final.fastq' \
  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
  --kit '3prime:v3' \
  --expected_cells 100 \
  --matrix_min_genes 10 \
  --matrix_min_cells 1 \
  --matrix_max_mito 50 \
  --out_dir 'result_nobc'


#Minimum without bc

nextflow run epi2me-labs/wf-single-cell \
  -r master \
  -profile singularity \
  --fastq 'nobcExtracted_Final.fastq' \
  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
  --kit '3prime:v3' \
  --expected_cells 1 \
  --matrix_min_genes 1 \
  --matrix_min_cells 1 \
  --matrix_max_mito 50 \
  --out_dir 'result_nobc_min'
