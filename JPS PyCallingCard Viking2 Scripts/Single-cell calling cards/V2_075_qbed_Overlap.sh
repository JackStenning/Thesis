#!/bin/bash

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/Results

module load BEDTools/2.30.0-GCC-11.2.0

#!/bin/bash

for dir in */; do
    cd "$dir" || continue
    for res in *result*/; do
        cd "$res" || continue
        for cc in *.ccf; do
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/Andy_q0.05/Andy_ChIP_q0.05_peaks.bed > "1kbp_A_${cc}_B_Andy.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4_q0.05/Encode_ER_ChIP_q0.05_peaks.bed > "1kbp_A_${cc}_B_BRD4.bed"
        done
        cd "$(dirname "$PWD")"
    done
    cd "$(dirname "$PWD")"
done

/mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/E1_barcodeFollowUpTest/initial_test_ccf_files

for res in *result*/; do
        cd "$res" || continue
        for cc in *.bed; do
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/Andy_q0.05/Andy_ChIP_q0.05_peaks.bed > "1kbp_A_${cc}_B_Andy.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4_q0.05/Encode_ER_ChIP_q0.05_peaks.bed > "1kbp_A_${cc}_B_BRD4.bed"
        done
	cd "$(dirname "$PWD")"
done