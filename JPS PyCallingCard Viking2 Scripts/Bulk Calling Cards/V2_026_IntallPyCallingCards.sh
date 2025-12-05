#!/bin/bash

module load Apptainer/latest
module load docker/6.0.1
module load Miniconda3/23.5.2-0

source activate CallingCardsPython

python3 -m pip install pyparsing twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy gnuplot
pip install "git+https://github.com/The-Mitra-Lab/pycallingcards.git" --upgrade
pip install git+https://github.com/cmatkhan/callingCardsTools.git