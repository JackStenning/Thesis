#!/bin/sh

srun --ntasks=1 --time=00:30:00 bedtools genomecov -ibam GSM4471638_HCT-116_HyPBase.bam -bg -scale 47 > GSM4471638_HCT-116_HyPBase.bedgraph

