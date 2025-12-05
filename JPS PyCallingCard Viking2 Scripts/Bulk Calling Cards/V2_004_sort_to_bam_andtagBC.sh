#!/usr/bin/env bash
#SBATCH --job-name=uptobed               # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-08:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=1               # Number of CPU cores per MPI task
#SBATCH --mem=100G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log


### bash sort to bam
#Purge then load modules
		module purge
		module load SAMtools/1.16.1-GCC-11.3.0
		
#Sorting to bam
cd 202306_JS_BulkCC/
for dir in *; do
  cd "$dir"
  echo "now in "$dir""
	for file in *.sam; do
	  echo "$file"
		samtools view \
 	  	 	-bS -h -F 260 \
 	  		 "$file" | \
 	   		samtools sort - -o "$file"_mapped.bam
		done
    cd ../
done
cd ../
### bash tag barcodes

#Purge then load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0

#Activate python environment
module load Miniconda3

cd 202306_JS_BulkCC/
command -v python
source activate CallingCardsPy
command -v python
python3 -m pip install twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy


#Tag barcodes
#Purge then load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0

#Activate python environment
module load Miniconda3/23.5.2-0
#Activate python environment

command -v python
source activate CallingCardsPy
command -v python
python3 -m pip install twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy

#Tag barcodes
cd 202306_JS_BulkCC/
for direct in ER*; do
cd "$direct"
	for file in *mapped.bam; do
		python ../../TagBam.py \
 		   --tag XP:Z:CGT \
 		   "$file" \
 		   "$file"_tagged.bam	
	done
	cd ../
done
for direct in WT*; do
cd "$direct"
	for file in *mapped.bam; do
		python ../../TagBam.py \
 		   --tag XP:Z:TAG \
 		   "$file" \
 		   "$file"_tagged.bam	
	done
	cd ../
done
cd ../
