import os
os.chdir('/mnt/scratch/users/jps558/Stenning_Data_Analysis/2024-JS-ssCC/CallingCards_ONT/E1_barcodeFollowUpTest/')

import pysam
from Bio.Seq import Seq
from difflib import SequenceMatcher
import regex

# === USER INPUTS ===

input_bam = "aligned.sorted.bam"

# Define adaptors
forward_adaptors = [
    "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AAAGATAGTCTGCGTAAAATTGACGC"
]

# Input/output
input_bam_path = "aligned.sorted.bam"
trimmed_bam_path = "trimmed_nobc.bam"
untrimmed_bam_path = "untrimmed_nobc.bam"
trimmed_fastq_path = "trimmed_nobc.fastq"
untrimmed_fastq_path = "untrimmed_nobc.fastq"

# Create list with forward + reverse complements
adaptors = forward_adaptors + [str(Seq(a).reverse_complement()) for a in forward_adaptors]

# Trimming behavior
MAX_MISMATCHES = 2        # 0 = exact match, 1 = strict fuzzy
EXTRACT_BARCODES = True
CB_LEN = 16               # Cell barcode length
UB_LEN = 12               # UMI length


# Compile regex patterns with fuzzy matching
compiled_adaptors = [regex.compile(f'({a}){{e<={MAX_MISMATCHES}}}') for a in adaptors]

# ============
# TRIM FUNCTION (fast + fuzzy match)
# ============

def trim_if_adaptor_found(read):
    seq = read.query_sequence
    qual = read.query_qualities

    if seq is None or qual is None:
        return read, False

    for pattern in compiled_adaptors:
        match = pattern.search(seq)
        if match:
            i = match.start()
            end = match.end()

            # Extract barcode/UMI if desired
            if EXTRACT_BARCODES:
                cb_start = end
                cb = seq[cb_start:cb_start + CB_LEN]
                ub = seq[cb_start + CB_LEN:cb_start + CB_LEN + UB_LEN]

                if len(cb) == CB_LEN:
                    read.set_tag("CB", cb)
                if len(ub) == UB_LEN:
                    read.set_tag("UB", ub)

            # Trim read (remove adaptor and anything after)
            trimmed_seq = seq[:i] + seq[end:]
            trimmed_qual = qual[:i] + qual[end:]

            read.query_sequence = trimmed_seq
            read.query_qualities = trimmed_qual
            return read, True

    return read, False

# ============
# FASTQ WRITER
# ============

def write_fastq(fh, read):
    fq_seq = read.query_sequence
    fq_qual = read.query_qualities
    if fq_seq is None or fq_qual is None:
        return
    fq_qual_str = "".join([chr(q + 33) for q in fq_qual])
    fh.write(f"@{read.query_name}\n{fq_seq}\n+\n{fq_qual_str}\n")

# ============
# MAIN SCRIPT
# ============

with pysam.AlignmentFile(input_bam_path, "rb") as in_bam, \
     pysam.AlignmentFile(trimmed_bam_path, "wb", template=in_bam) as trimmed_bam, \
     pysam.AlignmentFile(untrimmed_bam_path, "wb", template=in_bam) as untrimmed_bam, \
     open(trimmed_fastq_path, "w") as trimmed_fq, \
     open(untrimmed_fastq_path, "w") as untrimmed_fq:

    trimmed_count = 0
    total_count = 0

    for read in in_bam.fetch(until_eof=True):
        total_count += 1
        trimmed_read, was_trimmed = trim_if_adaptor_found(read)

        if was_trimmed:
            trimmed_bam.write(trimmed_read)
            write_fastq(trimmed_fq, trimmed_read)
            trimmed_count += 1
        else:
            untrimmed_bam.write(read)
            write_fastq(untrimmed_fq, read)

    print(f"Total reads processed: {total_count}")
    print(f"Reads trimmed: {trimmed_count}")
    print(f"Reads untrimmed: {total_count - trimmed_count}")