#!/bin/bash

cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT

###TEST

## Illumina primers

grep --ignore-case ACACTCTTTCCCTACACGACGCTCTTCCGATCT E1_barcode22/E1_barcode22_cat.fastq > fullIllumina_primer.txt
grep --ignore-case ACACGACGCTCTTCCGATCT E1_barcode22/E1_barcode22_cat.fastq > Illumina_primer_afterbioT.txt
grep --ignore-case AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT E1_barcode22/E1_barcode22_cat.fastq > fullIllumina_primer_revComp.txt
grep --ignore-case AGATCGGAAGAGCGTCGTGT E1_barcode22/E1_barcode22_cat.fastq > Illumina_primer_afterbioT.txt

### PB_LTR primer #Minus CTAG

grep --ignore-case gcgtcaattttacgcagactatcttt E1_barcode22/E1_barcode22_cat.fastq > fullLTR_primer.txt
grep --ignore-case aaagatagtctgcgtaaaattgacgc E1_barcode22/E1_barcode22_cat.fastq > fullLTR_primer_revComp.txt



#Degenerate barcode sequence
grep --ignore-case -E 'AGACTATCTTT[ATGC]{4}GGTTAA' E1_barcode22/E1_barcode22_cat.fastq > E1_barcodes.txt


grep --ignore-case -E 'TTAACC[ATGC]{4}AAAGATAGTCT' E1_barcode22/E1_barcode22_cat.fastq > E1_barcodes_revComp.txt

grep --ignore-case -E 'AGACTATCTTT[ATGC]{4}GGTTAA' fullIllumina_primer.txt > E1_barcodes_AfterPrimer.txt


grep --ignore-case -E 'TTAACC[ATGC]{4}AAAGATAGTCT' fullIllumina_primer_revComp.txt > E1_barcodes_AfterPrimer_RevComp.txt


#Scan all files

for dir in *barcode*/; do
	cd "$dir"
	echo "$dir"
	for reads in *.fastq; do
	echo "$reads"
	grep --ignore-case ACACTCTTTCCCTACACGACGCTCTTCCGATCT "$reads" > "$reads"_fullIllumina_primer.txt
	grep --ignore-case AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "$reads" > "$reads"_fullIllumina_primer_revComp.txt
	grep --ignore-case gcgtcaattttacgcagactatcttt "$reads" > "$reads"_fullLTR_primer.txt
	grep --ignore-case aaagatagtctgcgtaaaattgacgc "$reads" > "$reads"_fullLTR_primer_revComp.txt
	grep --ignore-case -E 'AGACTATCTTT[ATGC]{4}GGTTAA' "$reads"_fullIllumina_primer.txt > "$reads"_fullIllumina_primer_PB_TR_Filter.txt
	grep --ignore-case -E 'TTAACC[ATGC]{4}AAAGATAGTCT' "$reads"_fullIllumina_primer.txt > "$reads"_fullIllumina_primer_PB_TR_revComp_Filter.txt
	grep --ignore-case -E 'AGACTATCTTT[ATGC]{4}GGTTAA' "$reads"_fullIllumina_primer_revComp.txt > "$reads"_fullIllumina_primer_revComp_PB_TR_Filter.txt
	grep --ignore-case -E 'TTAACC[ATGC]{4}AAAGATAGTCT' "$reads"_fullIllumina_primer_revComp.txt > "$reads"_fullIllumina_primer_revComp_PB_TR_revComp_Filter.txt
	grep --ignore-case ACACTCTTTCCCTACACGACGCTCTTCCGATCT "$reads"_fullLTR_primer.txt > fullLTR_primer_And_Illumina.txt
	grep --ignore-case AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "$reads"_fullLTR_primer.txt > fullLTR_primer_revComp_And_Illumina.txt
	grep --ignore-case ACACTCTTTCCCTACACGACGCTCTTCCGATCT "$reads"_fullLTR_primer_revComp.txt > fullLTR_primer_revComp_And_Illumina.txt
	grep --ignore-case AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "$reads"_fullLTR_primer_revComp.txt > fullLTR_primer_revComp_revComp_And_Illumina.txt
	cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
	done
done


#Scan for full primer
cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
#Scan for LTR primer
for dir in *barcode*/; do
	cd "$dir"
	echo "$dir"
	for reads in *.fastq; do
	echo "$reads"
	grep --ignore-case GCGTCAATTTTACGCAGACTATCTTTCTAG "$reads" > "$reads"MismatchedLTR_primer.txt
	grep --ignore-case CTAGAAAGATAGTCTGCGTAAAATTGACGC "$reads" > "$reads"MismatchedLTR_primer_Rev_Comp.txt
	cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
	done
done

#Find reverse primer
for dir in *barcode*/; do
	cd "$dir"
	echo "$dir"
	for reads in *.fastq; do
	echo "$reads"
	grep --ignore-case AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "$reads"MismatchedLTR_primer.txt > "$reads"MismatchedLTR_primer_PlusIlluminaRevComp.txt
	grep --ignore-case ACACTCTTTCCCTACACGACGCTCTTCCGATCT "$reads"MismatchedLTR_primer_Rev_Comp.txt > "$reads"MismatchedLTR_primer_Rev_Comp_PlusIllumina.txt
	cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
	done
done
#Find primer LTR plus GGTTAA for full LTR 
for dir in *barcode*/; do
	cd "$dir"
	echo "$dir"
	for reads in *.fastq; do
	echo "$reads"
	grep --ignore-case GCGTCAATTTTACGCAGACTATCTTTCTAGggttaa "$reads" > "$reads"MismatchedLTR_primer_WithGGTTAA.txt
	grep --ignore-case ttaaccCTAGAAAGATAGTCTGCGTAAAATTGACGC "$reads" > "$reads"MismatchedLTR_primer_Rev_Comp_WithGGTTAA.txt
	cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
	done
done


#Find imperfect primer
for dir in *barcode*/; do
	cd "$dir"
	echo "$dir"
	for reads in *.fastq; do
	echo "$reads"
	grep -E "AAAGATAGTCTGCGTAAAATTGACGC.{0,3}$" "$reads" | grep -E "CACGACGCTCTTCCGATCT" > "$reads"AndyImperfectMatch.txt
	cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT
	done
done

grep -E "AAAGATAGTCTGCGTAAAATTGACGC.{0,3}$" | grep -E "CACGACGCTCTTCCGATCT"



