#!/bin/sh

srun --ntasks=1 --time=00:30:00 samtools sort SRR2481799_aligned.bam -o SRR2481799_sorted.bam
