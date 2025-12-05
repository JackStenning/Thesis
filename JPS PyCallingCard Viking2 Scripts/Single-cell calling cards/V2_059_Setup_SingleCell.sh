#!/bin/bash

cp -R  /mnt/scratch/projects/biol-holding-2023/2024-JS-ssCC/ /mnt/scratch/users/jps558/Stenning_Data_Analysis/

cd 2024-JS-ssCC/GeneExpression/

for dir in JS*/; do
	cd "$dir"
	echo "$dir"
	for file in *.fastq.gz; do
	gunzip "$file"
	done
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/GeneExpression/
done

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT

cd E1_barcode22/
cat *fastq.gz > E1_barcode22_cat.fastq.gz
gunzip E1_barcode22_cat.fastq.gz

cd ../E2_barcode23/
cat *fastq.gz > E2_barcode23_cat.fastq.gz
gunzip E2_barcode23_cat.fastq.gz

cd ../E3_barcode24/
cat *fastq.gz > E3_barcode24_cat.fastq.gz
gunzip E3_barcode24_cat.fastq.gz

cd ../W1_barcode18/
cat *fastq.gz > W1_barcode18_cat.fastq.gz
gunzip W1_barcode18_cat.fastq.gz

cd ../W2_barcode19/
cat *fastq.gz > W2_barcode19_cat.fastq.gz
gunzip W2_barcode19_cat.fastq.gz

cd ../W3_barcode20/
cat *fastq.gz > W3_barcode20_cat.fastq.gz
gunzip W3_barcode20_cat.fastq.gz

cd ../W4_barcode21/
cat *fastq.gz > W4_barcode21_cat.fastq.gz
gunzip W4_barcode21_cat.fastq.gz


### Run after quota exceeded
cd JS_W3/
for file in *.fastq.gz; do
gunzip "$file"
done

cd ../JS_W4/
for file in *.fastq.gz; do
gunzip "$file"
done

