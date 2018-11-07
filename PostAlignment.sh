
#!/bin/bash

#variables
File_Dir=/root/fastq_dir
Picard_Path=/github.com/broadinstitute/picard/releases/download/2.18.16
Bedtools_Path=/bedtools2/bin
INPUT=HISAT2_virome_alignment_first_pass_Cufflinks

for f in $(find $File_Dir/ -name *${INPUT}.sam); do
  FILENAME=$(echo $f | cut -d'_' -f1,2)
  #convert .sam to .bam
  samtools view -Su ${FILENAME}_${INPUT}.sam > ${FILENAME}_${INPUT}.bam

  #sort .bam
  samtools sort ${FILENAME}_${INPUT}.bam > ${FILENAME}_${INPUT}_sorted.bam

  #MarkDuplicates with Picard

  java -jar $Picard_Path/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${FILENAME}_${INPUT}_sorted.bam \
  O=${FILENAME}_${INPUT}_sorted_MrkDup.bam M=${FILENAME}_marked_dup_metrics.txt 

  #capture unique reads
  samtools view -bq 50 ${FILENAME}_${INPUT}_sorted_MrkDup.bam > ${FILENAME}_${INPUT}_sorted_MrkDup_unique.bam

  #remove tandem repeats
  ${Bedtools_Path}/bedtools intersect -abam ${FILENAME}_${INPUT}_sorted_MrkDup_unique.bam -b $File_Dir/Shortened_Virome_2.7.7.80.10.50.2000.bed -v > ${FILENAME}_${INPUT}_sorted_MrkDup_unique_noTanRep.bam

  #index
  samtools index ${FILENAME}_${INPUT}_sorted_MrkDup_unique_noTanRep.bam
  samtools index ${FILENAME}_${INPUT}_sorted_MrkDup_unique.bam
  samtools index ${FILENAME}_${INPUT}_sorted_MrkDup.bam
  samtools index ${FILENAME}_${INPUT}_sorted.bam

  #count reads
  samtools idxstats ${FILENAME}_${INPUT}_sorted_MrkDup_unique_noTanRep.bam > ${FILENAME}virome_counts_MrkDup_unique_noTanRep.txt
  samtools idxstats ${FILENAME}_${INPUT}_sorted_MrkDup_unique.bam > ${FILENAME}virome_counts_MrkDup_unique.txt
  samtools idxstats ${FILENAME}_${INPUT}_sorted_MrkDup.bam > ${FILENAME}virome_counts_MrkDup.txt
  samtools idxstats ${FILENAME}_${INPUT}_sorted.bam > ${FILENAME}virome_counts.txt
done