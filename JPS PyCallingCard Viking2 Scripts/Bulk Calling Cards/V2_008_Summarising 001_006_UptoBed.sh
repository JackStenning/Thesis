#!/bin/bash

### bash cut barcode
#Load modules
module load cutadapt/3.4-GCCcore-10.3.0


#First step is to cut off the barcode and surrounding adapter
#For ER calling cards in this experiment, the barcode is the same - CGT

cd 202306_JS_BulkCC/
for direct in ER*; do
cd "$direct"
gunzip *.gz
	for file in *1.fq; do
		cutadapt \
		    -g ^CGTTTTACGCAGACTATCTTTCTAGGGTTAA \
		    --minimum-length 1 \
		    --discard-untrimmed \
		    -e 0 \
		    --no-indels \
	 	   -o "$file"TrimmedBC.fq \
	 	   "$file"
	done
	cd ../
done
cd ../

#For WT calling cards in this experiment, the barcode is the same - TAG

cd 202306_JS_BulkCC/
for direct in WT*; do
cd "$direct"
gunzip *.gz
	for file in *1.fq; do
		cutadapt \
		    -g ^TAGTTTACGCAGACTATCTTTCTAGGGTTAA \
		    --minimum-length 1 \
		    --discard-untrimmed \
		    -e 0 \
		    --no-indels \
	 	   -o "$file"TrimmedBC.fq \
	 	   "$file"
	done
	cd ../
done
cd ../




###bash cut index
#Load modules
module load cutadapt/3.4-GCCcore-10.3.0

#Cut adapt will take degenerate sequence for indexes - NNNNNNNNNN
cd 202306_JS_BulkCC/
for dir in *; do
  cd "$dir"
  echo "Now in "$dir""
  for file in *TrimmedBC.fq; do
  	cutadapt \
  	-a CTGTCTCTTATACACATCTCCGAGCCCACGAGACTNNNNNNNNNNTCTCGTATGCCGTCTTCTGCTTG \
  	--minimum-length 1 \
  	-o "$file"_fulltrimmed.fastq \
  	"$file"
  done
  cd ../
done
cd ../

### bash align 

#purge then load modules
module purge
module load deepTools/3.5.0-foss-2021a
module load BWA/0.7.17-foss-2019b
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0

cd 202306_JS_BulkCC/
for dir in *; do
  cd "$dir"
  echo "Now in "$dir""
  for file in *fulltrimmed.fastq; do
   	echo "Processing "$file""
	  bwa mem -p -t 40 /mnt/scratch/users/jps558/staging/hg38/hg38.fa "$file" > "$file".sam
	  done
	  cd ../
done
cd ../




### bash sort to bam
#Purge then load modules
		module purge
		module load SAMtools/1.16.1-GCC-11.3.0
		
#Sorting to bam
cd 202306_JS_BulkCC/
for dir in *; do
  cd "$dir"
  echo "now in "$dir""
	for file in *.sam; do
	  echo "$file"
		samtools view \
 	  	 	-bS -h -F 260 \
 	  		 "$file" | \
 	   		samtools sort - -o "$file"_mapped.bam
		done
    cd ../
done
cd ../


### Set up env
cd 202306_JS_BulkCC/
python3 -m venv .CallingCardsPy
source .CallingCardsPy/bin/activate
python3 -m pip install pyparsing twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy
cd ../



### bash tag barcodes

#Purge then load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0

#Tag barcodes
#Purge then load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0

#Activate python environment
module load Miniconda3/23.5.2-0
#Activate python environment

command -v python
source activate CallingCardsPy
command -v python
python3 -m pip install twobitreader pysam numpy pandas scipy statsmodels pybedtools astropy

#Tag barcodes
cd 202306_JS_BulkCC/
for direct in ER*; do
cd "$direct"
	for file in *mapped.bam; do
		python ../../TagBam.py \
 		   --tag XP:Z:CGT \
 		   "$file" \
 		   "$file"_tagged.bam	
	done
	cd ../
done
for direct in WT*; do
cd "$direct"
	for file in *mapped.bam; do
		python ../../TagBam.py \
 		   --tag XP:Z:TAG \
 		   "$file" \
 		   "$file"_tagged.bam	
	done
	cd ../
done
cd ../

#TAG WT barcodes
echo "bash TAG WT barcodes"
cd 202306_JS_BulkCC/
for direct in WT*; do
cd "$direct"
	for file in *mapped.bam; do
		python ../../TagBam.py \
 		   --tag XP:Z:TAG \
 		   "$file" \
 		   "$file"_tagged.bam	
	done
	cd ../
done
cd ../


### Bash tag index
echo "bash tag index"
cd 202306_JS_BulkCC/
#ER1
cd ER1/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:ACTCGCTA \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../
#ER2
cd ER2/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:GGAGCTAC \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../
#ER3
cd ER3/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:GCGTAGTA \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../
#ER5
cd ER5/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:CGGAGCCT \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../
#ER6
cd ER6/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:TACGCTGC \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../
################################################################################################
#WT1
cd WT1/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:ATGCGCAG \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../
#WT2
cd WT2/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:TAGCGCTC \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../
#WT3
cd WT3/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:ACTGAGCG \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../
#WT4
cd WT4/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:CCTAAGAC \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../
#WT5
cd WT5/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:CGATCAGT \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../
#WT6
cd WT6/
for file in *_tagged.bam; do
	python ../../TagBam.py \
	   --tag XJ:Z:TGCAGCTA \
	   "$file" \
	   "$file"_tagged2.bam	
done
cd ../

### bash annotate insert
echo "bash annotate insert"
#Purge then load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0

#Prepare python environment
module load Miniconda3/23.5.2-0
source activate CallingcardsPy

cd 202306_JS_BulkCC
for dir in *; do
	cd "$dir"
	echo "$dir"
	for file in *tagged2.bam; do
		python ../../AnnotateInsertionSites.py \
		    --transposase PB \
		    -f \
		    "$file" \
		    /mnt/scratch/users/jps558/staging/hg38/hg38.2bit \
 		   "$file"_final.bam
		done
	for mile in *mapped.bam_tagged.bam_tagged2.bam_final.bam; do
            mv "$mile" "${mile%%mapped.bam_tagged.bam_tagged2.bam_final.bam}final.bam"
        done
	cd  ../
done
cd ../

### bash index bam
echo "bash index bam"
module purge
module load SAMtools/1.16.1-GCC-11.3.0

cd 202306_JS_BulkCC/
for dir in *; do
	cd "$dir"
	echo "$dir"
	for file in *final.bam; do
		samtools index "$file"
		done
	cd  ../
done
cd ../

### bash convert bam to CC
echo "bash convert bam to CC"
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0

module load Miniconda3/23.5.2-0
source activate CallingcardsPy

cd 202306_JS_BulkCC/
mkdir outputs
mkdir outputs/ER
mkdir outputs/WT
for dir in ER*; do
	cd "$dir"
	for file in *final.bam; do
		python ../../BamToCallingCard.py \
		    -b XP XJ \
		    -i "$file" \
		    > "$file".qbed
		done
	for mile in *final.bam.qbed; do
            mv "$mile" "${mile%%final.bam.qbed}.qbed"
        done
  for save in *.qbed; do
    cp *.qbed ../outputs/ER
    done
  for store in *final.bam*; do
    cp "$store" ../outputs/ER
    done
	cd  ../
done

for dir in WT*; do
	cd "$dir"
	for file in *final.bam; do
		python ../../BamToCallingCard.py \
		    -b XP XJ \
		    -i "$file" \
		    > "$file".qbed
		done
	for mile in *final.bam.qbed; do
            mv "$mile" "${mile%%final.bam.qbed}.qbed"
        done
  for save in *.qbed; do
    cp *.qbed ../outputs/WT
    done
  for store in *final.bam*; do
    cp "$store" ../outputs/WT
    done
	cd  ../
done

###TO FIX CONCANTENATED FILES, I'LL BE STARTING FROM CONVERTING INTO QBEDS

cd 202306_JS_BulkCC/outputs/ER
for file in *.qbed; do
	mv "$file" "${file%%.fqTrimmedBC.fq_fulltrimmed.fastq.sam_.qbedfinal.qbed}.qbed"
	done
for dile in *final.bam.bai; do
	mv "$dile" "${dile%%.fqTrimmedBC.fq_fulltrimmed.fastq.sam_final.bam.baifinal.bam.bai}_final.bam.bai"
	done
for mile in *final.bam; do
	mv "$mile" "${mile%%.fqTrimmedBC.fq_fulltrimmed.fastq.sam_final.bamfinal.bam}_final.bam"
	done


module load SAMtools/1.16.1-GCC-11.3.0
AND BED TOOLS


#Cat ER
mkdir combined/
cat ER*HF552DSX7_L2_1_.qbed | bedtools sort -i > combined/ER_HyPB_HF552DSX7_L2_1.qbed
cat ER*HFK3FDSX7_L2_1_.qbed | bedtools sort -i > combined/ER_HyPB_HFK3FDSX7_L2_1.qbed
cat ER*.qbed | bedtools sort -i > combined/ER_HyPB_EKDL230008858.qbed
samtools merge ER_HyPB_EKDL230008858_final.bam ER*final.bam
cd combined/
mv *qbed ../

#Cat WT
mkdir combined/
cat WT*HF552DSX7_L2_1_.qbed | bedtools sort -i > combined/WT_HyPB_HF552DSX7_L2_1.qbed
cat WT*HFK3FDSX7_L2_1_.qbed | bedtools sort -i > combined/WT_HyPB_HFK3FDSX7_L2_1.qbed
cat WT*.qbed | bedtools sort -i > combined/WT_HyPB_EKDL230008858.qbed
samtools merge WT_HyPB_EKDL230008858_final.bam WT*final.bam
cd combined/
mv *qbed ../
cd ../../../../
