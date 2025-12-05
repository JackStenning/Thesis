#!/usr/bin/env bash
#SBATCH --job-name=Renser               # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=0-10:00:00               # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=1                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=1               # Number of CPU cores per MPI task
#SBATCH --mem=200G                        # Total memory to request
#SBATCH --account=biol-cards-2023        # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jps558@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log


module load R
Rscript -e "rmarkdown::render('Viking SP1 and ESR1 bulk calling card analysis.Rmd')"
