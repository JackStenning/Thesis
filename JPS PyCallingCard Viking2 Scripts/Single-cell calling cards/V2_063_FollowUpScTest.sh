cd /mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/FollowUpTest

###TEST

## Illumina primers

grep --ignore-case ACACTCTTTCCCTACACGACGCTCTTCCGATCT JS_ES1_allpass_adapter_trimmed.fastq > fullIllumina_primer.txt
grep --ignore-case AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT JS_ES1_allpass_adapter_trimmed.fastq > fullIllumina_primer_revComp.txt

grep --ignore-case aaagatagtctgcgtaaaattgacgc fullIllumina_primer.txt > fullIllumina_primer_withLTR_Rev.txt

grep --ignore-case gcgtcaattttacgcagactatcttt fullIllumina_primer_revComp.txt > fullIllumina_primer_revComp_withLTR.txt




### PB_LTR primer #Minus CTAG

grep --ignore-case gcgtcaattttacgcagactatcttt JS_ES1_allpass_adapter_trimmed.fastq > fullLTR_primer.txt
grep --ignore-case aaagatagtctgcgtaaaattgacgc JS_ES1_allpass_adapter_trimmed.fastq > fullLTR_primer_revComp.txt

#Andy Imperfect match
grep -E "AAAGATAGTCTGCGTAAAATTGACGC.{0,3}$" JS_ES1_allpass_adapter_trimmed.fastq | grep -E "CACGACGCTCTTCCGATCT" > JS_ES1_allpass_adapter_trimmed.fastq_AndyImperfectMatch.txt

#TTAA match

grep -E "TTAACC[ATGC]{4}AAAGATAGTCTGCGTAAAATTGACGC.{0,3}$" JS_ES1_allpass_adapter_trimmed.fastq | grep -E "CACGACGCTCTTCCGATCT" > JS_ES1_allpass_adapter_trimmed.fastqTTAAImperfectMatch.txt

#TTAA full match

grep -E "TTAACC[ATGC]{4}AAAGATAGTCTGCGTAAAATTGACGC.{0,3}$" JS_ES1_allpass_adapter_trimmed.fastq | grep -E "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" > JS_ES1_allpass_adapter_trimmed.fastq_PB_TTAARevComp_IlluminaFwd_.txt
	
grep -E "gcgtcaattttacgcagactatcttt[ATGC]{4}GGTTAA.{0,3}$" JS_ES1_allpass_adapter_trimmed.fastq | grep -E "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" > JS_ES1_allpass_adapter_trimmed.fastq_PB_TTAAFwd_IlluminaRevComp_.txt

#TTAA full match - no anchor

grep -E "TTAACC[ATGC]{4}AAAGATAGTCTGCGTAAAATTGACGC" JS_ES1_allpass_adapter_trimmed.fastq | grep -E "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" > Non_Ending_PB_TTAARevComp_IlluminaFwd_.txt
	
grep -E "gcgtcaattttacgcagactatcttt[ATGC]{4}GGTTAA" JS_ES1_allpass_adapter_trimmed.fastq | grep -E "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" > Non_Ending_PB_TTAAFwd_IlluminaRevComp_.txt

