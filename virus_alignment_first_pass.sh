#!/bin/bash

#variables
File_Dir=/root/fastq_dir

# run!
for f in $(find $File_Dir/ -name "*forward_paired.fq"); do
  FILENAME=$(echo $f | cut -d'_' -f1,2)
  hisat2 --phred33 --known-splicesite-infile $File_Dir/Shortened_Virome_ss.txt --dta-cufflinks \
     -x $File_Dir/genome_Virome_index -1 ${FILENAME}_forward_paired.fq -2 ${FILENAME}_reverse_paired.fq \
     -U ${FILENAME}_forward_unpaired.fq,${FILENAME}_reverse_unpaired.fq \
     -S ${FILENAME}_HISAT2_virome_alignment_first_pass_Cufflinks.sam \
     --un ${FILENAME}_virome_first_pass_Cufflinks_un-seqs --un-conc ${FILENAME}_virome_first_pass_Cufflinks_un-conc-mate
done