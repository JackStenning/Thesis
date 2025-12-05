#!/bin/bash

module load gnuplot/5.4.8-GCCcore-12.3.0
module load SAMtools/1.17-GCC-12.2.0

mkdir /mnt/scratch/users/jps558/Stenning_Data_Analysis/CC_Pipeline/bamstats
mkdir /mnt/scratch/users/jps558/Stenning_Data_Analysis/CC_Pipeline/bamstats/test1
mkdir /mnt/scratch/users/jps558/Stenning_Data_Analysis/CC_Pipeline/bamstats/test2

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/CC_Pipeline/Pipelineoutput/ER_AllReps/hops/samtools

plot-bamstats -p /mnt/scratch/users/jps558/Stenning_Data_Analysis/CC_Pipeline/bamstats/test1/test1 ER_AllReps_failing_merged_sorted.stats

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/CC_Pipeline/Pipelineoutput2/ER_AllReps/hops/samtools

plot-bamstats -p /mnt/scratch/users/jps558/Stenning_Data_Analysis/CC_Pipeline/bamstats/test2/test2 ER_AllReps_failing_merged_sorted.stats
