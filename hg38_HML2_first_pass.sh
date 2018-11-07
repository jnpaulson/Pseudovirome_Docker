#!/bin/bash

#variables
File_Dir=/root/fastq_dir
OUTPUT=HISAT2_hg38_HML2_alignment

for f in $(find $File_Dir/ -name "*forward_paired.fq"); do
  FILENAME=$(echo $f | cut -d'_' -f1,2)
  hisat2 --phred33 --known-splicesite-infile ${File_Dir}/hg38_HML2_genes_051918_ss.txt --dta-cufflinks \
    -x ${File_Dir}/genome_hg38_HML2_index -1 ${FILENAME}_forward_paired.fq -2 ${FILENAME}_reverse_paired.fq \
    -U ${FILENAME}_forward_unpaired.fq,${FILENAME}_reverse_unpaired.fq -S ${FILENAME}_${OUTPUT}.sam --un ${FILENAME}_${OUTPUT}_un-seqs \ 
    --un-conc ${FILENAME}_${OUTPUT}_un-conc-mate;
  cat ${FILENAME}_${OUTPUT}_un-conc-mate.1 ${FILENAME}_${OUTPUT}_un-conc-mate.2 ${FILENAME}_${OUTPUT}_un-seqs > ${FILENAME}_${OUTPUT}_combined.fq
done