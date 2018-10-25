#!/bin/bash
#SBATCH --job-name=hg38_back_alignment   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH -c 1                  # Run a single task
#SBATCH --mem-per-cpu=10gb           # Memory per node
#SBATCH --time=6-23:59:59             # Time limit hrs:min:sec
#SBATCH --output=hg38_round2_%A-%a.out    # Standard output and error log
#SBATCH --array=0-106                 # Array range

#load modules
module purge
module load centos6/0.0.1-fasrc01
module load slurm-drmaa/1.2.0-fasrc01
module load bedtools2/2.26.0-fasrc01
module load samtools/1.5-fasrc02
module load hisat2/2.1.0-fasrc01

#list all .fq files in directory and print their names
FILE=($(ls *HISAT2_virome_alignment_first_pass_Cufflinks.sam |rev | cut -c 49- | rev | uniq))
echo ${FILE}

# grab out filename from the array exported from our 'parent' shell
FILENAME=${FILE[$SLURM_ARRAY_TASK_ID]}
echo ${FILENAME}

#sam to bam
samtools view -Su ${FILENAME}HISAT2_virome_alignment_first_pass_Cufflinks.sam > ${FILENAME}HISAT2_virome_alignment_first_pass_Cufflinks.bam

#sort bam
samtools sort ${FILENAME}HISAT2_virome_alignment_first_pass_Cufflinks.bam > ${FILENAME}HISAT2_virome_alignment_first_pass_Cufflinks_sorted.bam

#bam to fastq
bedtools bamtofastq -i ${FILENAME}HISAT2_virome_alignment_first_pass_Cufflinks_sorted.bam -fq ${FILENAME}_R1.fq -fq2 ${FILENAME}_R2.fq

#align pulled reads to hg38
hisat2 --phred33 --known-splicesite-infile /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/froy/hg38_with_HML2_genome_files/hg38_HML2_genes_051918_ss.txt --dta-cufflinks  -x /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/froy/hg38_with_HML2_genome_files/genome_hg38_HML2_index -1 ${FILENAME}_R1.fq -2 ${FILENAME}_R2.fq -S ${FILENAME}_HISAT2_hg38_HML2_alignment_second_pass.sam --un ${FILENAME}_HISAT2_hg38_HML2_alignment_second_pass_un-seqs --un-conc ${FILENAME}_HISAT2_hg38_HML2_alignment_second_pass_un-conc-mate


