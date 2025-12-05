#!/usr/bin/env bash
#SBATCH --job-name=Nextflowpipe          # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-20:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=12               # Number of CPU cores per MPI task
#SBATCH --mem=256G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log


cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/Combined

### Load Modules ###

#CutAdapt Modules
module load cutadapt/4.2
module load SeqKit/2.3.1
module load seqtk/1.3-GCC-11.3.0

#Epi2Me/Nextflow modules
module load Nextflow/23.04.2
module load Apptainer/latest

#Load nextflow pipeline
nextflow run epi2me-labs/wf-artic --help

#Update workflow
nextflow pull epi2me-labs/wf-single-cell -r master


#Calling Card Modules
module load BEDTools/2.31.0-GCC-12.3.0
module load BEDOPS/2.4.41-foss-2021b
module load SAMtools/1.16.1-GCC-11.3.0
module load Apptainer/latest
module load Python/3.11.3-GCCcore-12.3.0

### Trim Primers ###

## Cut Adapt
# Define primers

FORWARD="ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
REVERSE_nobc="AAAGATAGTCTGCGTAAAATTGACGC"

Combine files

cat E1_barcode22_cat.fastq.cleancombined.fastq E2_barcode23_cat.fastq.cleancombined.fastq E3_barcode24_cat.fastq.cleancombined.fastq > Ecat.trimmed.cleancombined.fastq

cat W1_barcode18_cat.fastq.cleancombined.fastq W2_barcode19_cat.fastq.cleancombined.fastq W1_barcode20_cat.fastq.cleancombined.fastq W4_barcode21_cat.fastq.cleancombined.fastq > Wcat.trimmed.cleancombined.fastq

# Start loop

for dir in *barcode*/; do
	cd $dir
	for fq in *cat.trimmed.cleancombined.fastq; do

		#Remove duplicates
		seqkit rmdup -s $fq > $fq.finalcombined.fastq

## Nextflow

		nextflow run epi2me-labs/wf-single-cell \
		  -r master \
		  -profile singularity \
		  --fastq "${fq}.finalcombined.fastq" \
		  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
		  --kit '3prime:v3' \
		  --expected_cells 100 \
		  --matrix_min_genes 10 \
		  --matrix_min_cells 1 \
		  --matrix_max_mito 50 \
		  --out_dir "${fq}result_trim_full"

		nextflow run epi2me-labs/wf-single-cell \
		  -r master \
		  -profile singularity \
		  --fastq "${fq}.finalcombined.fastq" \
		  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
		  --kit '3prime:v3' \
		  --expected_cells 1 \
		  --matrix_min_genes 1 \
		  --matrix_min_cells 1 \
		  --matrix_max_mito 50 \
		  --out_dir "${fq}result_trim_minimum"
		done
	for res in *result*/; do
		cd $res
		cd *cat*/
		for bam in *cat*.tagged.bam; do
			python /mnt/scratch/users/jps558/staging/AnnotateInsertionSites.py \
				--transposase PB \
				"${bam}" \
				/mnt/scratch/users/jps558/staging/hg38/hg38.2bit \
				"${bam}.CC_taged.bam"

			cp "${bam}.bai" "${bam}.CC_taged.bam.bai" 

			python /mnt/scratch/users/jps558/staging/UMIFilter.py \
			    -p 10x \
			    -i "${bam}.CC_taged.bam" \
			    --verbose \
			    -o "${bam}.filtered.final.bam"

			cp "${bam}.CC_taged.bam.bai" "${bam}.filtered.final.bam.bai"

			python /mnt/scratch/users/jps558/staging/BamToCallingCard.py \
			    -b CB \
			    -i "${bam}.CC_taged.bam" \
			    > "${bam}.unsorted.ccf"

			grep -v '^chrM' "${bam}.unsorted.ccf" > "${bam}_final.ccf"
			sort-bed "${bam}_final.ccf" > "${bam}_final.sorted.ccf"

			python /mnt/scratch/users/jps558/staging/BamToCallingCard.py \
			    -b CB \
			    -i "${bam}.filtered.final.bam" \
			    > "${bam}_UMIFilt.unsorted.ccf"
			
			sort-bed "${bam}_UMIFilt.unsorted.ccf" > "${bam}_UMIFilt_final.ccf"
			
			mkdir /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/$dir
			mkdir /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/$dir$res
			
			cp "${bam}_UMIFilt_final.ccf" /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/$dir$res
			cp "${bam}_final.sorted.ccf" /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/$dir$res
			done
		cd "$(dirname "$PWD")"
		cd "$(dirname "$PWD")"
		done
	/mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/Combined
done




