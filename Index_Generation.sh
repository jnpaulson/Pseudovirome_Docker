#!/bin/bash
Ref_Dir=/root/fastqc_dir
#generate .ss and .exon files
/hisat2-2.1.0/extract_splice_sites.py $Ref_Dir/Shortened_Virome.gtf > $Ref_Dir/Shortened_Virome_ss.txt
/hisat2-2.1.0/extract_exons.py $Ref_Dir/Shortened_Virome.gtf > $Ref_Dir/Shortened_Virome_exons.txt
/hisat2-2.1.0/extract_splice_sites.py $Ref_Dir/hg38_HML2_genes_051918.gtf > $Ref_Dir/hg38_HML2_genes_051918_ss.txt
/hisat2-2.1.0/extract_exons.py $Ref_Dir/hg38_HML2_genes_051918.gtf > $Ref_Dir/hg38_HML2_genes_051918_exons.txt


#generate  indexes
/hisat2-2.1.0/hisat2-build --ss $Ref_Dir/Shortened_Virome_ss.txt --exon $Ref_Dir/Shortened_Virome_exons.txt -f $Ref_Dir/Shortened_Virome.fa $Ref_Dir/genome_Virome_index
/hisat2-2.1.0/hisat2-build --ss $Ref_Dir/hg38_HML2_genes_051918_ss.txt --exon $Ref_Dir/hg38_HML2_genes_051918_exons.txt -f $Ref_Dir/genome.fa  $Ref_Dir/genome_hg38_HML2_index