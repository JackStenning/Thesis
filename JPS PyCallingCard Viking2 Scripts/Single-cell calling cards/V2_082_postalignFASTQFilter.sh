#!/usr/bin/env bash
#SBATCH --job-name=SetupAltPipefq          # Job name
#SBATCH --partition=nodes                  # Partition
#SBATCH --time=0-25:00:00                  # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                         # Number of MPI tasks
#SBATCH --cpus-per-task=12                 # Number of CPU cores per task
#SBATCH --mem=256G                         # Total memory
#SBATCH --account=biol-cards-2023          # Project account
#SBATCH --mail-type=END,FAIL               # Mail events
#SBATCH --mail-user=jps558@york.ac.uk      # Email
#SBATCH --output=%x-%j.log                 # Standard output log
#SBATCH --error=%x-%j.err                  # Standard error log

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/E1_barcodeFollowUpTest/

#Load modules
module purge
module load GCCcore/12.3.0
module load Python/3.10.8-GCCcore-12.2.0
module load BEDTools/2.31.0-GCC-12.3.0
module load minimap2/2.26-GCCcore-12.3.0
module load SAMtools/1.20-GCC-13.2.0
module load BEDOPS/2.4.41-foss-2021b
#Create environment 
python -m venv ~/myenv
source ~/myenv/bin/activate

Upgrade pip and install pure-Python dependencies
pip install --upgrade pip
pip install pysam biopython regex twobitreader numpy pandas scipy statsmodels astropy pybedtools

#Confirm environment works
python -c 'import pybedtools, twobitreader, pysam; print("Environment OK!")'

# Process each FASTQ file
for file in trimmed*bc.fastq; do
    base="${file%.fastq}"

    # Convert FASTQ to SAM
    minimap2 -ax map-ont /mnt/scratch/users/jps558/staging/hg38/hg38.fa "${file}" > "${base}.sam"

    # Convert SAM to BAM, sort, and index
    samtools view -Sb "${base}.sam" > "${base}_fq.bam"
    samtools sort "${base}_fq.bam" -o "${base}_fq.sorted.bam"
    samtools index "${base}_fq.sorted.bam"

    # Annotate insertion sites
    python /mnt/scratch/users/jps558/staging/AnnotateInsertionSites.py \
        --transposase PB \
        "${base}_fq.sorted.bam" \
        /mnt/scratch/users/jps558/staging/hg38/hg38.2bit \
        "${base}_fq.CC_tagged.bam"

    cp "${base}_fq.bai" "${base}_fq.CC_tagged.bam.bai"

    # Filter UMIs
    python /mnt/scratch/users/jps558/staging/UMIFilter.py \
        -p 10x \
        -i "${base}_fq.CC_tagged.bam" \
        --verbose \
        -o "${base}_fq.filtered.final.bam"

    cp "${base}_fq.CC_tagged.bam.bai" "${base}_fq.filtered.final.bam.bai"

    # Generate CallingCard files
    python /mnt/scratch/users/jps558/staging/BamToCallingCard.py \
        -b CB \
        -i "${base}_fq.CC_tagged.bam" \
        > "${base}_fq.unsorted.ccf"

    grep -v '^chrM' "${base}_fq.unsorted.ccf" > "${base}_fq_final.ccf"
    sort-bed "${base}_fq_final.ccf" > "${base}_fq_final.sorted.ccf"

    python /mnt/scratch/users/jps558/staging/BamToCallingCard.py \
        -b CB \
        -i "${base}_fq.filtered.final.bam" \
        > "${base}_fq_UMIFilt.unsorted.ccf"

    sort-bed "${base}_fq_UMIFilt.unsorted.ccf" > "${base}_fq_UMIFilt_final.ccf"

    # Copy results
    cp "${base}_fq_UMIFilt_final.ccf" PostAlignOutput/
    cp "${base}_fq_final.sorted.ccf" PostAlignOutput/
done
