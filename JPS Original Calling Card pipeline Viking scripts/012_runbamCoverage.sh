#!/bin/sh

srun --ntasks=1 --time=00:30:00 bamCoverage -b GSM4471638_HCT-116_HyPBase.bam -of bigwig -o GSM4471638_HCT-116_HyPBase.bw
