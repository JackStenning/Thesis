#!/bin/sh

srun --ntasks=1 --time=00:30:00 samtools view -S -b SRR2481799_aligned.sam > SRR2481799_aligned.bam

