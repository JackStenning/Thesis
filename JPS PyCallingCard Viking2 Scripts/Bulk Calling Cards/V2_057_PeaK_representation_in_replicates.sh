#!/bin/bash
srun --time=03:00:00 --partition=interactive --mem=10GB --pty /bin/bash

module load BEDTools/2.30.0-GCC-11.2.0

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis

cp -R Replicate_analysis /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/Replicate_analysis/combined/

for rep in filtered*.qbed; do
#	bedtools window -w 1000 -c -a /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_2/peak_data_ER_WT2.bed -b $rep > ../rep_overlap/a_peakset2_b_$rep.bed
#	bedtools window -w 1000 -c -a /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_4/peak_data_ER_WT4.bed -b $rep > ../rep_overlap/a_peakset4_b_$rep.bed
	bedtools window -u -a $rep -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_4/peak_data_ER_WT4.bed > ../rep_overlap/a_"$rep"_b_peakset4_.bed
	bedtools window -u -a $rep -b /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/peakset_2/peak_data_ER_WT2.bed > ../rep_overlap/a_"$rep"_b_peakset2_.bed

done
