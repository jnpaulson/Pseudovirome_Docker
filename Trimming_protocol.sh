
#!/bin/bash

# grab out filename from list and trim

Trimming_Path=/Trimmomatic-0.36/adapters
FASTQC_Path=/FastQC

for f in $(find /fastq -name "*R1_001.fastq"); do
  FILENAME=$(echo $f | cut -d'_' -f1,2)
  $FASTQC_Path/fastqc ${FILENAME}_R1_001.fastq ${FILENAME}_R2_001.fastq 
  java -jar /Trimmomatic-0.36/trimmomatic-0.36.jar  PE \
  -phred33 ${FILENAME}_R1_001.fastq ${FILENAME}_R2_001.fastq ${FILENAME}_forward_paired.fq ${FILENAME}_forward_unpaired.fq ${FILENAME}_reverse_paired.fq ${FILENAME}_reverse_unpaired.fq \
  ILLUMINACLIP:$Trimming_Path/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:75 
  $FASTQC_Path/fastqc ${FILENAME}_forward_paired.fq ${FILENAME}_forward_unpaired.fq ${FILENAME}_reverse_paired.fq ${FILENAME}_reverse_unpaired.fq
done