#!/bin/bash

module load BEDTools/2.30.0-GCC-11.2.0

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_2025Analysis/

#Define variables
result=/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/results/
root=/mnt/scratch/users/jps558/Stenning_Data_Analysis/Final_year/2024_25Analysis/

#Change directory to pvalue files
for dir in *ChIP/; do
	cd "$dir"
	for mac in MACS2*/; do
		cd "$mac"
		for pval in */; do
			cd "$pval"
			cp "$root"*peak*


