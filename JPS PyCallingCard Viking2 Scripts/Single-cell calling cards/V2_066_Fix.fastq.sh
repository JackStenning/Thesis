#!/bin/bash

cp -R /mnt/scratch/projects/biol-holding-2023/2024-JS-ssCC/CallingCards_ONT/E1_barcode22 /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/E1_barcode22

cat *.fastq.gz > E1_barcode22_cat.fastq.gz

gunzip E1_barcode22_cat.fastq.gz

