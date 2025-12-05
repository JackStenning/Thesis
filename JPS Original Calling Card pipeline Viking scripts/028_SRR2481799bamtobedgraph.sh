#!/bin/sh

srun --ntasks=1 --time=00:30:00 bedtools genomecov -ibam SRR2481799_sorted.bam -bg -scale 35 > SRR2481799_sorted.bedgraph
