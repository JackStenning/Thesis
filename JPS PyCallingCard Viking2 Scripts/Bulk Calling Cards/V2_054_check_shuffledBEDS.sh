#!/bin/bash


############ ALL FILES BELOW CREATE OVERLAP OF 4-2%, AS SUCH A RANDOM FILE WAS SELECTED - a RANDOM NUMBER BETWEEN 1 AND 5 WAS GIVEN - THIS WAS 2, FROM ASCENDING ORDER THIS MEANS 155780119 WAS SELECTED ########



module load BEDTools/2.30.0-GCC-11.2.0

#Directory already has peak files

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/ROC

## Andy

bedtools shuffle -i Andy_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -chrom -noOverlapping -seed 862157398 > Arandom_sites_862157398.bed

bedtools shuffle -i Andy_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -chrom -noOverlapping -seed 932536857 > Arandom_sites_932536857.bed

bedtools shuffle -i Andy_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -chrom -noOverlapping -seed 463449987 > Arandom_sites_463449987.bed

bedtools shuffle -i Andy_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -chrom -noOverlapping -seed 33386894 > Arandom_sites_33386894.bed

bedtools shuffle -i Andy_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -chrom -noOverlapping -seed 155780119 > Arandom_sites_155780119.bed

## Encode

bedtools shuffle -i Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -chrom -noOverlapping -seed 862157398 > Erandom_sites_862157398.bed

bedtools shuffle -i Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -chrom -noOverlapping -seed 932536857 > Erandom_sites_932536857.bed

bedtools shuffle -i Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -chrom -noOverlapping -seed 463449987 > Erandom_sites_463449987.bed

bedtools shuffle -i Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -chrom -noOverlapping -seed 33386894 > Erandom_sites_33386894.bed

bedtools shuffle -i Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -chrom -noOverlapping -seed 155780119 > Erandom_sites_155780119.bed



cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/ROC

#Extend peak size by 1000 bp to avoid using window
bedtools slop -i peak_data_ER_WT2.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -b 1000 > slop_ER_WT2.bed

bedtools slop -i peak_data_ER_WT4.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -b 1000 > slop_ER_WT4.bed

bedtools slop -i Andy_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -b 1000 > slop_Andy_ChIP_q0.05_peaksNP.bed

bedtools slop -i Encode_ER_rep2_ChIP_q0.05_peaksNP.bed -g /mnt/scratch/users/jps558/staging/hg38/hg38.chrom.sizes -b 1000 > slop_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed


mkdir Arando

for filt in slop*.bed; do
	for ran in Arandom_sites*.bed; do
	bedtools intersect -a "$filt" -b "$ran" -wa -u > Arando/A_"$filt"_B_"$ran".bed
done
done

mkdir Erando

for filt in slop*.bed; do
	for ran in Erandom_sites*.bed; do
	bedtools intersect -a "$filt" -b "$ran" -wa -u > Erando/A_"$filt"_B_"$ran".bed
done
done

