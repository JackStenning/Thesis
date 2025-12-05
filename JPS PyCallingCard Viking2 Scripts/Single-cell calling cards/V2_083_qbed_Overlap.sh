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
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peak_data_ER_WT2.bed > "1kbp_A_${cc}_B_highCC.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peak_data_ER_WT4.bed > "1kbp_A_${cc}_B_lowCC.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38/bedpe_output/Chiapet_hg38_combined.bed_anchor1.bed > "1kbp_A_${cc}_B_chiaanchor1.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38/bedpe_output/Chiapet_hg38_combined.bed_anchor2.bed > "1kbp_A_${cc}_B_chiaanchor2.bed"
        done
        for cc in Full*data_ER_test2.bed; do
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/Andy_q0.05/Andy_ChIP_q0.05_peaks.bed > "1kbp_A_${cc}_B_Andy.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4_q0.05/Encode_ER_ChIP_q0.05_peaks.bed > "1kbp_A_${cc}_B_BRD4.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peak_data_ER_WT2.bed > "1kbp_A_${cc}_B_highCC.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peak_data_ER_WT4.bed > "1kbp_A_${cc}_B_lowCC.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38/bedpe_output/Chiapet_hg38_combined.bed_anchor1.bed > "1kbp_A_${cc}_B_chiaanchor1.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38/bedpe_output/Chiapet_hg38_combined.bed_anchor2.bed > "1kbp_A_${cc}_B_chiaanchor2.bed"
        done
        for cc in Full*data_ER_test4.bed; do
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/Andy_q0.05/Andy_ChIP_q0.05_peaks.bed > "1kbp_A_${cc}_B_Andy.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4_q0.05/Encode_ER_ChIP_q0.05_peaks.bed > "1kbp_A_${cc}_B_BRD4.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peak_data_ER_WT2.bed > "1kbp_A_${cc}_B_highCC.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peak_data_ER_WT4.bed > "1kbp_A_${cc}_B_lowCC.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38/bedpe_output/Chiapet_hg38_combined.bed_anchor1.bed > "1kbp_A_${cc}_B_chiaanchor1.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38/bedpe_output/Chiapet_hg38_combined.bed_anchor2.bed > "1kbp_A_${cc}_B_chiaanchor2.bed"
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
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peak_data_ER_WT2.bed > "1kbp_A_${cc}_B_highCC.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peak_data_ER_WT4.bed > "1kbp_A_${cc}_B_lowCC.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38/bedpe_output/Chiapet_hg38_combined.bed_anchor1.bed > "1kbp_A_${cc}_B_chiaanchor1.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38/bedpe_output/Chiapet_hg38_combined.bed_anchor2.bed > "1kbp_A_${cc}_B_chiaanchor2.bed"
        done
	cd "$(dirname "$PWD")"
done

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/E1_barcodeFollowUpTest/FinalCCF/Post\ align\ before\ fix

        for cc in *.bed; do
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/Andy_q0.05/Andy_ChIP_q0.05_peaks.bed > "1kbp_A_${cc}_B_Andy.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/BRD4_q0.05/Encode_ER_ChIP_q0.05_peaks.bed > "1kbp_A_${cc}_B_BRD4.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peak_data_ER_WT2.bed > "1kbp_A_${cc}_B_highCC.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peak_data_ER_WT4.bed > "1kbp_A_${cc}_B_lowCC.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38/bedpe_output/Chiapet_hg38_combined.bed_anchor1.bed > "1kbp_A_${cc}_B_chiaanchor1.bed"
            bedtools window -w 1000 -c -a "$cc" -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/chiapet/chia_pet_hg38/bedpe_output/Chiapet_hg38_combined.bed_anchor2.bed > "1kbp_A_${cc}_B_chiaanchor1.bed"
        done
