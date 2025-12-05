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


cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT

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

# Start loop

for dir in *FollowUpTest/; do
	cd $dir
	for fq in JS_ES1_allpass_adapter_trimmed.fastq; do
		## This is trying to make the reads more like true scCC reads before trimming to help downstream steps.
		#Cut out anything that doesnt have compelte primer match without PB barcode

		cutadapt \
		  -g "$FORWARD;min_overlap=33" \
		  -a "$REVERSE_nobc;min_overlap=26" \
		  --no-indels \
		  --rc \
		  -e 0 \
		  --minimum-length 1 \
		  --debug \
		  --discard-untrimmed \
		  -o $fq.trimmed.fastq \
		  $fq

		# Reverse complimenting (to get PB TR at 5') and trimming 
		seqtk seq -r $fq.trimmed.fastq > RC$fq.trimmed.fastq

		#Clip TTAA
		cutadapt \
		 -g ^NNNNGGTTAA \
		 --no-indels \
		 -e 0 \
		 --minimum-length 1 \
		 --discard-untrimmed \
		 -o $fq.trimmed1.fastq \
		 RC$fq.trimmed.fastq

		#Trim not RC incase of missing reads
		cutadapt \
		 -g ^NNNNGGTTAA \
		 --no-indels \
		 -e 0 \
		 --minimum-length 1 \
		 --discard-untrimmed \
		 -o $fq.trimmed2.fastq \
		 $fq.trimmed.fastq

		#Combine
		cat $fq.trimmed1.fastq $fq.trimmed2.fastq > $fq.trimmedcombined.fastq

		#Remove duplicates
		seqkit rmdup -s $fq.trimmedcombined.fastq > $fq.cleancombined.fastq

		cp $fq.cleancombined.fastq /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results/Combined
		
## Nextflow

		nextflow run epi2me-labs/wf-single-cell \
		  -r master \
		  -profile singularity \
		  --fastq "${fq}.cleancombined.fastq" \
		  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
		  --kit '3prime:v3' \
		  --expected_cells 100 \
		  --matrix_min_genes 10 \
		  --matrix_min_cells 1 \
		  --matrix_max_mito 50 \
		  --out_dir 'result_trim_full'

		nextflow run epi2me-labs/wf-single-cell \
		  -r master \
		  -profile singularity \
		  --fastq "${fq}.cleancombined.fastq" \
		  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
		  --kit '3prime:v3' \
		  --expected_cells 1 \
		  --matrix_min_genes 1 \
		  --matrix_min_cells 1 \
		  --matrix_max_mito 50 \
		  --out_dir 'result_trim_minimum'

		nextflow run epi2me-labs/wf-single-cell \
		  -r master \
		  -profile singularity \
		  --fastq "${fq}" \
		  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
		  --kit '3prime:v3' \
		  --expected_cells 100 \
		  --matrix_min_genes 10 \
		  --matrix_min_cells 1 \
		  --matrix_max_mito 50 \
		  --out_dir 'result_full'

		nextflow run epi2me-labs/wf-single-cell \
		  -r master \
		  -profile singularity \
		  --fastq "${fq}" \
		  --ref_genome_dir '/mnt/scratch/users/jps558/staging/hg38/' \
		  --kit '3prime:v3' \
		  --expected_cells 1 \
		  --matrix_min_genes 1 \
		  --matrix_min_cells 1 \
		  --matrix_max_mito 50 \
		  --out_dir 'result_minimum'
		done
	for res in result*/; do
		cd $res
		cd *allpass*/
		for bam in *allpass*trimmed.tagged.bam; do
			python /mnt/scratch/users/jps558/staging/AnnotateInsertionSites.py \
				--transposase PB \
				"${bam}" \
				/mnt/scratch/users/jps558/staging/hg38/hg38.2bit \
				"${bam}.CC_tagged.bam"

			cp "${bam}.bai" "${bam}.CC_tagged.bam.bai" 

			python /mnt/scratch/users/jps558/staging/UMIFilter.py \
			    -p 10x \
			    -i "${bam}.CC_tagged.bam" \
			    --verbose \
			    -o "${bam}.filtered.final.bam"

			cp "${bam}.CC_tagged.bam.bai" "${bam}.filtered.final.bam.bai"

			python /mnt/scratch/users/jps558/staging/BamToCallingCard.py \
			    -b CB \
			    -i "${bam}.CC_tagged.bam" \
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
	cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
done



