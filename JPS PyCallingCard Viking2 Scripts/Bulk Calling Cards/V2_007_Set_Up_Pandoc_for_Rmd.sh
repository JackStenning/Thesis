#!/bin/sh
cd 202306_JS_BulkCC
wget https://github.com/jgm/pandoc/releases/download/3.1.9/pandoc-3.1.9-linux-amd64.tar.gz
tar -xzf pandoc-3.1.9-linux-amd64.tar.gz
cd pandoc-3.1.9/
cd bin
mkdir ~/bin
mv pandoc ~/bin

#TO RUN
module load R
#Rscript -e "rmarkdown::render('Viking SP1 and ESR1 bulk calling card analysis.Rmd')"

cd ../