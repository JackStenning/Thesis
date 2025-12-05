#!/usr/bin/env bash
#SBATCH --job-name=cutadapt               # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-08:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=1               # Number of CPU cores per MPI task
#SBATCH --mem=50G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log


### bash cut barcode
#Load modules
module load cutadapt/3.4-GCCcore-10.3.0


#First step is to cut off the barcode and surrounding adapter
#For ER calling cards in this experiment, the barcode is the same - CGT

cd 202306_JS_BulkCC/
for direct in ER*; do
cd "$direct"
gunzip *.gz
	for file in *1.fq; do
		cutadapt \
		    -g ^CGTTTTACGCAGACTATCTTTCTAGGGTTAA \
		    --minimum-length 1 \
		    --discard-untrimmed \
		    -e 0 \
		    --no-indels \
	 	   -o "$file"TrimmedBC.fq \
	 	   "$file"
	done
	cd ../
done
cd ../

#For WT calling cards in this experiment, the barcode is the same - TAG

cd 202306_JS_BulkCC/
for direct in WT*; do
cd "$direct"
gunzip *.gz
	for file in *1.fq; do
		cutadapt \
		    -g ^TAGTTTACGCAGACTATCTTTCTAGGGTTAA \
		    --minimum-length 1 \
		    --discard-untrimmed \
		    -e 0 \
		    --no-indels \
	 	   -o "$file"TrimmedBC.fq \
	 	   "$file"
	done
	cd ../
done
cd ../




###bash cut index
#Load modules
module load cutadapt/3.4-GCCcore-10.3.0

#Cut adapt will take degenerate sequence for indexes - NNNNNNNNNN
cd 202306_JS_BulkCC/
for dir in *; do
  cd "$dir"
  echo "Now in "$dir""
  for file in *TrimmedBC.fq; do
  	cutadapt \
  	-a CTGTCTCTTATACACATCTCCGAGCCCACGAGACTNNNNNNNNNNTCTCGTATGCCGTCTTCTGCTTG \
  	--minimum-length 1 \
  	-o "$file"_fulltrimmed.fastq \
  	"$file"
  done
  cd ../
done
cd ../

sbatch V2_002_Alignment.sh
