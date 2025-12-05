#!/bin/bash

#Goto Directory
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/FollowUpTest

##### TEST #####

#Load module
module load cutadapt/4.2
module load SeqKit/2.3.1
module load seqtk/1.3-GCC-11.3.0



# Define primers
FORWARD="ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
REVERSE_nobc="AAAGATAGTCTGCGTAAAATTGACGC"
REVERSE_bc="TTAACCNNNNAAAGATAGTCTGCGTAAAATTGACGC"

# Input and output
INPUT="JS_ES1_allpass_adapter_trimmed.fastq"
OUTPUT_1="extracted_nobc.fastq"
OUTPUT_2="extracted_bc.fastq"
ROUGH="rough_extracted.fastq"

# Use cutadapt to find and extract all matching internal regions


#Below is 17/07

cutadapt \
  -g "$FORWARD...$REVERSE_bc;min_overlap=33" \
  --no-indels \
  --rc \
  -e 0 \
  --debug \
  --discard-untrimmed \
  -o Extracted.fastq \
  "$INPUT"

cutadapt \
  -g "$FORWARD...$REVERSE_bc;min_overlap=33" \
  --no-indels \
  --rc \
  -e 0 \
  --debug \
  --minimum-length 1 \
  --discard-untrimmed \
  -o Extracted_min1.fastq \
  "$INPUT"


#############################################
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
  -o nobcExtracted.fastq \
  "$INPUT"

# Reverse complimenting (to get PB TR at 5') and trimming 
seqtk seq -r nobcExtracted.fastq > RC_nobcExtracted.fastq

#Clip TTAA
cutadapt \
 -g ^NNNNGGTTAA \
 --no-indels \
 -e 0 \
 --minimum-length 1 \
 --discard-untrimmed \
 -o RC_nobcExtracted_TTAA.fastq \
 RC_nobcExtracted.fastq

#Trim not RC incase of missing reads
cutadapt \
 -g ^NNNNGGTTAA \
 --no-indels \
 -e 0 \
 --minimum-length 1 \
 --discard-untrimmed \
 -o nobcExtracted_TTAA.fastq \
 nobcExtracted.fastq

#Combine
cat RC_nobcExtracted_TTAA.fastq nobcExtracted_TTAA.fastq > nobcExtracted_TTAAcombined.fastq

#Remove duplicates
seqkit rmdup -s nobcExtracted_TTAAcombined.fastq > nobcExtracted_Final.fastq


