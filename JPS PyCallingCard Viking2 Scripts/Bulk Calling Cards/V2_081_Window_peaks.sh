#!/bin/bash

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year

module load BEDTools/2.30.0-GCC-11.2.0

# Only keep first 3 columns (chrom, start, end)
awk 'BEGIN{OFS="\t"} {print $1, int($2), int($3)}' BRD4/MACS2/Control_q0.05/BRD4_control_q0.05_peaksNP.bed > BRD4/MACS2/Control_q0.05/BRD4_control_q0.05peaks.bed
awk 'BEGIN{OFS="\t"} {print $1, int($2), int($3)}' BRD4_Zheng/MACS2/Control_Z_q0.05/BRD4_control_Z_q0.05_peaksNP.bed > BRD4_Zheng/MACS2/Control_Z_q0.05/BRD4_control_Z_q0.05peaks.bed

##Define comparisons
#BRD4
BRD4=BRD4/MACS2/Control_q0.05/BRD4_control_q0.05peaks.bed

# BRD4Z
BRD4Z=BRD4_Zheng/MACS2/Control_Z_q0.05/BRD4_control_Z_q0.05peaks.bed

mkdir overlap
mkdir overlap/cc_BRD4
mkdir overlap/cc_Chia
mkdir overlap/Chia_cc

#Calling Card overlap

for file in *eak*.bed; do
bedtools window -w 1000 -c -a "${file}" -b "$BRD4" > overlap/cc_BRD4/count_BRD4_"${file}"
bedtools window -w 1000 -c -a "${file}" -b "$BRD4Z" > overlap/cc_BRD4/count_BRD4Z_"${file}"
    for chia in chiapet/chia_pet_hg38/bedpe_output/*anchor*; do
        base=$(basename "$file" .bed)    # strip extension
        chia_base=$(basename "$chia" .bed)
        bedtools window -w 1000 -c -a "$file" -b "$chia" > overlap/cc_Chia/count_Chia_${base}_${chia_base}.txt
        bedtools window -w 1000 -a "$file" -b "$chia" > overlap/cc_Chia/basic_Chia_${base}_${chia_base}.txt
        bedtools window -w 1000 -u -a "$file" -b "$chia" > overlap/cc_Chia/uniq_basic_Chia_${base}_${chia_base}.txt
    done
    for chia in chiapet/chia_pet_hg38/Chiapet*.bed; do
    base=$(basename "$file" .bed)    # strip extension
    chia_base=$(basename "$chia" .bed)
    bedtools window -w 1000 -a "$chia" -b "$file" > overlap/Chia_cc/basic_Chia_${chia_base}_${base}.txt
    bedtools window -w 1000 -a "$chia" -b "$file" -u > overlap/Chia_cc/uniq_basic_Chia_${chia_base}_${base}.txt
    done
done

#ChIP overlap
mkdir overlap/chip_BRD4
mkdir overlap/chip_Chia

bedtools window -w 1000 -c -a Andy_q0.05/Andy_ChIP_q0.05_peaksNP.bed -b "$BRD4" > overlap/chip_BRD4/count_BRD4_Andy_ChIP_q0.05_peaksNP.bed
bedtools window -w 1000 -c -a Andy_q0.05/Andy_ChIP_q0.05_peaksNP.bed -b "$BRD4Z" > overlap/chip_BRD4/count_BRD4Z_Andy_ChIP_q0.05_peaksNP.bed
bedtools window -w 1000 -c -a BRD4_q0.05/Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -b "$BRD4" > overlap/chip_BRD4/count_BRD4_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed
bedtools window -w 1000 -c -a BRD4_q0.05/Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -b "$BRD4Z" > overlap/chip_BRD4/count_BRD4Z_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed
for chia in chiapet/chia_pet_hg38/bedpe_output/*anchor*; do
        chia_base=$(basename "$chia" .bed)
	bedtools window -w 1000 -c -a Andy_q0.05/Andy_ChIP_q0.05_peaksNP.bed -b "$chia" > overlap/chip_Chia/count_Chia_Andy_ChIP_q0.05_peaksNP.bed${chia_base}.txt
	bedtools window -w 1000 -a Andy_q0.05/Andy_ChIP_q0.05_peaksNP.bed -b "$chia" > overlap/chip_Chia/basic_Chia_Andy_ChIP_q0.05_peaksNP.bed${chia_base}.txt
	bedtools window -w 1000 -u -a Andy_q0.05/Andy_ChIP_q0.05_peaksNP.bed -b "$chia" > overlap/chip_Chia/uniq_basic_Chia_Andy_ChIP_q0.05_peaksNP.bed${chia_base}.txt
	bedtools window -w 1000 -c -a BRD4_q0.05/Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -b "$chia" > overlap/chip_Chia/count_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed${chia_base}.txt
	bedtools window -w 1000 -a BRD4_q0.05/Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -b "$chia" > overlap/chip_Chia/basic_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed${chia_base}.txt
	bedtools window -w 1000 -u -a BRD4_q0.05/Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -b "$chia" > overlap/chip_Chia/uniq_basic_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed${chia_base}.txt
done

#Check BEDPE
BEDPE_DIR="chiapet/chia_pet_hg38/bedpe_output"
OUTPUT_DIR="overlap/BEDPE"

mkdir -p "$OUTPUT_DIR"

# Loop over all BEDPE {file}s
BEDPE_DIR="chiapet/chia_pet_hg38/bedpe_output"
OUTPUT_DIR="overlap/BEDPE"

mkdir -p "$OUTPUT_DIR"

for FILE in *eak*.bed; do
    BASENAME_FILE=$(basename "$FILE" .bed)

    # Loop over all BEDPE files
    for BEDPE in "$BEDPE_DIR"/*.bedpe; do
        BASENAME_BEDPE=$(basename "$BEDPE" .bedpe)
        echo "Processing $BASENAME_BEDPE for $BASENAME_FILE"

        # Generate anchors
        ANCHOR1="$OUTPUT_DIR/${BASENAME_BEDPE}_anchor1.bed"
        ANCHOR2="$OUTPUT_DIR/${BASENAME_BEDPE}_anchor2.bed"

        awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' "$BEDPE" > "$ANCHOR1"
        awk 'BEGIN{OFS="\t"}{print $4,$5,$6}' "$BEDPE" > "$ANCHOR2"

        # List of datasets to check overlaps against
        DATASETS=(
            "$FILE"
            "Andy_q0.05/Andy_ChIP_q0.05_peaksNP.bed"
            "BRD4_q0.05/Encode_ER_rep2_ChIP_q0.05_peaksNP.bed"
        )

        # Loop through each dataset and run bedtools window
        for DATASET in "${DATASETS[@]}"; do
            BASENAME_DATASET=$(basename "$DATASET" .bed)
            
            bedtools window -w 1000 -a "$DATASET" -b "$ANCHOR1" \
                > "$OUTPUT_DIR/${BASENAME_BEDPE}_anchor1_${BASENAME_DATASET}_overlap.bed"
            
            bedtools window -w 1000 -a "$DATASET" -b "$ANCHOR2" \
                > "$OUTPUT_DIR/${BASENAME_BEDPE}_anchor2_${BASENAME_DATASET}_overlap.bed"
        done
    done
done

echo "All overlaps completed. Results in $OUTPUT_DIR"
