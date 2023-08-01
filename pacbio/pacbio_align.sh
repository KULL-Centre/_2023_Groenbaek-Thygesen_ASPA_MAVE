#!/bin/bash

# NOTE: samtools view -d np:X and -D np:X did not with version 1.11 but works with samtools-1.16.1
# stat the number of complete PacBio CCS passes of the insert
samtools view m54329U_201001_062023.ccs.bam | awk '{split($14,a,":"); print a[3]}' | sort -n | uniq -c > np_stat.txt

# Filter CCS reads to have 10 or more passes (tag: np)
seq 10 1000 > np.txt
samtools view -D np:np.txt  m54329U_201001_062023.ccs.bam | samtools fastq -T np > m54329U_201001_062023.np10.fastq
gzip m54329U_201001_062023.np10.fastq

# Make index for BWA alignment
bwa index aspa_ref_bwa.fasta

# Align
# bwa mem -t 4 aspa_ref_bwa.fasta m54329U_201001_062023.Q20.fastq.gz  1> aligned.sam  2> bwa.out
bwa mem -t 4 aspa_ref_bwa.fasta m54329U_201001_062023.np10.fastq.gz  1> aligned.sam  2> bwa.out
wc -l aligned.sam

# Filter out unaligned (4) and supplementary alignments (2048)
samtools view -F 2052 aligned.sam > aligned_filtered.sam
wc -l aligned_filtered.sam

# Make a fasta of aligned sequences
# Don't use samtools fasta, this will read FALG 16 and generate reverse-complement reads
awk '{if ($3 == "barcode-gfp-aspa") {print ">" $1 "  " $6 "\n" $10}}' aligned_filtered.sam > aligned_filtered.fasta

# Find barcodes (18 bases)
cutadapt -j 4 -e 0.2 -a ACAAATAGTT...TGCGAGTAGT -o barcodes.fasta aligned_filtered.fasta > cutadapt_bar.out

# Find variants (939 + 3 stop codons = 948 bases)
cutadapt -j 4 -e 0.2 -a AGCCACCATG...CTTAAGAATT -o variants.fasta aligned_filtered.fasta > cutadapt_var.out

# Back to one line per read
../tools/fasta2seq.sh barcodes.fasta > barcodes.seq
../tools/fasta2seq.sh variants.fasta > variants.seq

wc -l barcodes.seq
wc -l variants.seq

# Merge barcodes and variants using white-space
paste -d "  " barcodes.seq variants.seq > bar_var_raw.seq

# Check that reads are the same and print those that have the right length
awk '{ if ($1 == $3 && length($2)==18 && length($4)==948) print $2 "  " substr($4,1,939) }' bar_var_raw.seq > bar_var_filtered.txt
wc -l bar_var_filtered.txt

# Count unique barcode-variant pairs
sort bar_var_filtered.txt | uniq -c > bar_var_unq.txt
wc -l bar_var_unq.txt
