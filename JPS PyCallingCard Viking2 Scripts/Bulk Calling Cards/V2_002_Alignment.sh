#!/usr/bin/env bash
#SBATCH --job-name=Alignment               # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-20:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=1               # Number of CPU cores per MPI task
#SBATCH --mem=150G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log

### bash align 

#purge then load modules
module purge
module load deepTools/3.5.0-foss-2021a
module load BWA/0.7.17-foss-2019b
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0

cd 202306_JS_BulkCC/
for dir in *; do
  cd "$dir"
  echo "Now in "$dir""
  for file in *fulltrimmed.fastq; do
   	echo "Processing "$file""
	  bwa mem -p -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa "$file" > "$file".sam
	  done
	  cd ../
done
cd ../



