#!/bin/bash

### Set up env
cd 202306_JS_BulkCC/
python3 -m venv .CallingCardsPy
source .CallingCardsPy/bin/activate
python3 -m pip install pyparsing twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy
cd ../
