#!/bin/bash

#variables
File_Dir=/root/fastq_dir
Bedtools_Path=/bedtools2/bin
INPUT=HISAT2_virome_alignment_first_pass_Cufflinks
OUTPUT=HISAT2_hg38_HML2_alignment_second_pass

for f in $(find $File_Dir/ -name "*${INPUT}.bam"); do
  FILENAME=$(echo $f | cut -d'_' -f1,2)
  #sort .bam by name
  samtools sort -n ${FILENAME}_${INPUT}.bam > ${FILENAME}_${INPUT}_sorted_names.bam
  #extracting mapped reads from pseudovirome alignment
  samtools view -h -F ${FILENAME}_${INPUT}_sorted_names.bam |  samtools view -Su -  > ${FILENAME}_${INPUT}_sorted_names_mappedReads.bam
  #bam to fastq
  ${Bedtools_Path}/bedtools bamtofastq -i ${FILENAME}_${INPUT}_sorted_names_mappedReads.bam -fq ${FILENAME}_R1.fq -fq2 ${FILENAME}_R2.fq
  #align pulled reads to hg38
  hisat2 --phred33 --known-splicesite-infile ${File_Dir}/hg38_HML2_genes_051918_ss.txt --dta-cufflinks  \
    -x ${File_Dir}/genome_hg38_HML2_index -1 ${FILENAME}_R1.fq -2 ${FILENAME}_R2.fq -S ${FILENAME}_${OUTPUT}.sam --un ${FILENAME}_${OUTPUT}_un-seqs --un-conc ${FILENAME}_${OUTPUT}_un-conc-mate
done