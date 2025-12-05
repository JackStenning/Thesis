#!/usr/bin/env bash
#SBATCH --job-name=BRD4Chip         # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-20:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=12               # Number of CPU cores per MPI task
#SBATCH --mem=150G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log

#load modules
module load SRA-Toolkit/3.2.0-gompi-2024a

mkdir /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year

##curl -L -o GSE18046_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE18046&format=file"

#tar -xvf GSE18046_RAW.tar

#gunzip *.gz

#download interact files
wget -r -np -nH --cut-dirs=4 ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeGisChiaPet/ -P chia_pet/ -A "*.bed","*.bed.gz"

cd chia_pet
rm *Pol2*
rm *Ctcf*
rm *K562*

gunzip *gz

#Files were lifted over using the following tool - https://genome.ucsc.edu/cgi-bin/hgLiftOver
#Renamed ChiaPetRep1,2,3

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year

# Function to convert UCSC ChIA-PET BED to BEDPE

# Directories
INPUT_DIR="/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38"
OUTPUT_DIR="$INPUT_DIR/bedpe_output"

mkdir -p "$OUTPUT_DIR"

# Step 1: Combine all lifted BED replicates into a single file
COMBINED_BED="$OUTPUT_DIR/Chiapet_hg38_combined.bed"
cat "$INPUT_DIR"/ChiapetRep*.bed > "$COMBINED_BED"
echo "Combined BED file created: $COMBINED_BED"

# Step 2: Process each replicate + combined file
for FILE in "$INPUT_DIR"/ChiapetRep*.bed "$COMBINED_BED"; do
    BASENAME=$(basename "$FILE" "_lifted.bed")
    BEDPE="$OUTPUT_DIR/${BASENAME}.bedpe"

    # Create standard 6-column BEDPE (line 2 pasted onto line 1)
    awk 'NR % 2 == 1 {line1=$1"\t"$2"\t"$3; next} NR % 2 == 0 {print line1"\t"$1"\t"$2"\t"$3}' "$FILE" > "$BEDPE"

    # Generate anchor files
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$BEDPE" > "$OUTPUT_DIR/${BASENAME}_anchor1.bed"
    awk 'BEGIN{OFS="\t"} {print $4,$5,$6}' "$BEDPE" > "$OUTPUT_DIR/${BASENAME}_anchor2.bed"

    echo "Processed $BASENAME"
done