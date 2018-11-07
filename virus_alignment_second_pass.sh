#!/bin/bash
#variables
File_Dir=/root/fastq_dir
INPUT=HISAT2_hg38_HML2_alignment
OUTPUT=HISAT2_virome_alignment_Cufflinks

# run!
for f in $(find $File_Dir/ -name "*forward_paired.fq"); do
  FILENAME=$(echo $f | cut -d'_' -f1,2)
  hisat2 --phred33 --known-splicesite-infile ${File_Dir}/Shortened_Virome_ss.txt --dta-cufflinks \ 
      -x ${File_Dir}/genome_Virome_index  -q ${FILENAME}_${INPUT}_combined.fq -S ${FILENAME}_${OUTPUT}.sam --un ${FILENAME}_${OUTPUT}_un-seqs \
      --un-conc ${FILENAME}_${OUTPUT}_un-conc-mate
done