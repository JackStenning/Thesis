#!/bin/bash

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/202306_JS_BulkCC/ChIP/Check
cp -R /mnt/scratch/projects/biol-cards-2023/ERbigwigs /mnt/scratch/users/jps558/Stenning_Data_Analysis/202306_JS_BulkCC/ChIP/Check
cp -R Kent/ ERbigwigs/
cd ERbigwigs/
./Kent/bigWigMerge ER_A_90_SRR6652020_1.bw ER_B_90_SRR6652023_1.bw ER_C_90_SRR6652026_1.bw ER_D_90_SRR6652029_1.bw ER_combined_90_SRR665202.bg
./Kent/bedGraphToBigWig ER_combined_90_SRR665202.bg  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes ER_combined_90_SRR665202.bw
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis