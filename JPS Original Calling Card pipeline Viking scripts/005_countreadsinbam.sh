#!/bin/sh

srun --ntasks=1 --time=00:30:00 samtools view -c GSM4471638_HCT-116_HyPBase.bam

printf unmapped reads

srun --ntasks=1 --time=00:30:00 samtools view -c -F 4 GSM4471638_HCT-116_HyPBase.bam

printf mapped reads
