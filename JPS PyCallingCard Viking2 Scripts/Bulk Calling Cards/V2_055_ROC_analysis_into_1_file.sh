#!/bin/bash
srun --time=03:00:00 --partition=interactive --mem=10GB --pty /bin/bash

module load BEDTools/2.30.0-GCC-11.2.0

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year


#cp -R ROC backupROC

cd ROC/

#Extend peak size by 1000 bp to avoid using window
bedtools slop -i peak_data_ER_WT2.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -b 1000 > slop_ER_WT2.bed

bedtools slop -i peak_data_ER_WT4.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -b 1000 > slop_ER_WT4.bed

#Prepare for ROC calculation
mkdir output

#Filter based on p-value
for peak in slop_ER_WT*.bed; do 

        # Filter peaks based on current threshold T
	awk '$15 >= 0' "$peak" > output/filtered_"$peak"_0.bed
	awk '$15 >= 1e-13' "$peak" > output/filtered_"$peak"_13.bed
	awk '$15 >= 1e-12' "$peak" > output/filtered_"$peak"_12.bed
	awk '$15 >= 1e-11' "$peak" > output/filtered_"$peak"_11.bed
	awk '$15 >= 1e-10' "$peak" > output/filtered_"$peak"_10.bed
	awk '$15 >= 1e-9' "$peak" > output/filtered_"$peak"_9.bed
	awk '$15 >= 1e-8' "$peak" > output/filtered_"$peak"_8.bed
	awk '$15 >= 1e-7' "$peak" > output/filtered_"$peak"_7.bed
	awk '$15 >= 1e-6' "$peak" > output/filtered_"$peak"_6.bed
	awk '$15 >= 1e-5' "$peak" > output/filtered_"$peak"_5.bed
	awk '$15 >= 1e-3' "$peak" > output/filtered_"$peak"_3.bed
	awk '$15 >= 1' "$peak" > output/filtered_"$peak"_1.bed

done

#Change to output file and prep it
cd output/

mkdir andy
mkdir encode

#find TP and FP
for filt in filtered_slop*.bed; do

	#Find True positive in Andy and Encode ChIP
	bedtools intersect -a "$filt" -b ../Andy_ChIP_q0.05_peaksNP.bed -wa -u > andy/"$filt"_true_positive_andy.bed
	bedtools intersect -a "$filt" -b ../Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -wa -u > encode/"$filt"_true_positive_encode.bed

	
	#Count the lines
	wc -l andy/"$filt"_true_positive_andy.bed >> andy/ALLcount_true_positive_andy.bed
	wc -l encode/"$filt"_true_positive_encode.bed >> encode/ALLcount_true_positive_encode.bed

	#Find false positive

	bedtools subtract -a "$filt" -b ../Andy_ChIP_q0.05_peaksNP.bed -A > andy/"$filt"_false_positive_andy_.bed
	bedtools subtract -a "$filt" -b ../Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -A > encode/"$filt"_false_positive_encode_.bed


	#Count the lines
	wc -l andy/"$filt"_false_positive_andy_.bed >> andy/ALLcount_false_positive_andy_.bed
	wc -l encode/"$filt"_false_positive_encode_.bed >> encode/ALLcount_false_positive_encode_.bed


	#Find False negatives 
	bedtools subtract -a ../Andy_ChIP_q0.05_peaksNP.bed -b andy/"$filt"_true_positive_andy.bed -A > andy/"$filt"_false_negative_andy_.bed
	bedtools subtract -a ../Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -b encode/"$filt"_false_positive_encode_.bed -A > encode/"$filt"_false_negative_encode_.bed
	
	#Count the lines
	wc -l andy/"$filt"_false_negative_andy_.bed >> andy/ALLcount_false_negative_andy_.bed
	wc -l encode/"$filt"_false_negative_encode_.bed >> encode/ALLcount_false_negative_encode_.bed

	#Find True Negatives
	bedtools subtract -a Arandom_sites_155780119.bed -b "$filt" -A > andy/"$filt"_true_negative_andy.bed
	bedtools subtract -a Erandom_sites_155780119.bed -b "$filt" -A > encode/"$filt"_true_negative_encode.bed

	#Count the lines
	wc -l andy/"$filt"_true_negative_andy.bed >> andy/ALLcount_true_negative_andy.bed
	wc -l encode/"$filt"_true_negative_encode.bed >> encode/ALLcount_true_negative_encode.bed
done



	