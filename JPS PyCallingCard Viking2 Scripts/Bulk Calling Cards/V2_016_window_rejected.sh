#!/bin/bash

cd Rejectedpeaks/
module load BEDTools/2.30.0-GCC-11.2.0
mkdir results
for CC in *reject*.bed; do
	for ChIP in *GSM29704*; do
		echo $CC
		echo $ChIP
		bedtools window -c -w 1000 -a $CC -b $ChIP > results/1kbp_A_"$CC"_B_"$ChIP"
		echo "Next!"
		done
	done
cd ../