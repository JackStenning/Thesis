#!/bin/bash
cd Rejectedpeaks/results
module load R/4.2.1-foss-2022a
Rscript --vanilla ../ESR1_rejected.R
cd ../