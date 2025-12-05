#!/bin/bash
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
            mv "$mile" "${mile%%.fqTrimmedBC.fq_fulltrimmed.fastq.sam_mapped.bam_tagged.bam_tagged2.bam_final.bam}final.bam"
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
    cp *.qbed ../outputs/ER
    done
  for store in *final.bam*; do
    cp "$store" ../outputs/ER
    done
	cd  ../
done

#Cat ER
cd outputs/ER
cat ER*HF552DSX7_L2_1_.qbed | bedtools sort -i > ER_HyPB_HF552DSX7_L2_1.qbed
cat ER*HFK3FDSX7_L2_1_.qbed | bedtools sort -i > ER_HyPB_HFK3FDSX7_L2_1.qbed
cat ER*.qbed | bedtools sort -i > ER_HyPB_EKDL230008858.qbed
samtools merge ER*final.bam ER_HyPB_EKDL230008858_final.bam

#Cat WT
cd ../outputs/WT
cat WT*HF552DSX7_L2_1_.qbed | bedtools sort -i > WT_HyPB_HF552DSX7_L2_1.qbed
cat WT*HFK3FDSX7_L2_1_.qbed | bedtools sort -i > WT_HyPB_HFK3FDSX7_L2_1.qbed
cat WT*.qbed | bedtools sort -i > WT_HyPB_EKDL230008858.qbed
samtools merge WT*final.bam WT_HyPB_EKDL230008858_final.bam
cd ../../../

### bash cleanup 


for dir in 202306_JS_BulkCC/*; do
	cd "$dir"
	for file in *Trimmed*; do
		rm "$file"
		done
	for dile in *.sam; do
		rm "$dile"
		done
	for mile in tagged*.bam; do
		rm "$mile"
		done
	cd ../../
	done

