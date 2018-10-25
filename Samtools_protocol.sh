#!/bin/bash
#SBATCH --job-name=Samtools_PA   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH -c 1                  # Run a single task
#SBATCH --mem-per-cpu=10gb           # Memory per node
#SBATCH --time=6-23:59:59             # Time limit hrs:min:sec
#SBATCH --output=virus_samtools_%A-%a.out    # Standard output and error log
#SBATCH --array=0-106%10                 # Array range

#load modules
module purge
module load centos6/0.0.1-fasrc01
module load slurm-drmaa/1.2.0-fasrc01
module load samtools/1.5-fasrc02

#list all .sam virome alignment files in directory
FILE=($(ls *HISAT2_virome_alignment_Cufflinks.sam  | rev | cut -c 38- | rev | uniq)) 
 
# grab out filename from the array exported from our 'parent' shell
FILENAME=${FILE[$SLURM_ARRAY_TASK_ID]}
echo ${FILENAME}

#convert .sam to .bam
samtools view -Su ${FILENAME}HISAT2_virome_alignment_Cufflinks.sam > ${FILENAME}HISAT2_virome_alignment_Cufflinks.bam

#sort .bam
samtools sort ${FILENAME}HISAT2_virome_alignment_Cufflinks.bam > ${FILENAME}HISAT2_virome_alignment_Cufflinks_sorted.bam

#grab MAPQ50 reads
samtools view -bq 50 ${FILENAME}HISAT2_virome_alignment_Cufflinks_sorted.bam > ${FILENAME}HISAT2_virome_alignment_Cufflinks_sorted_unique.bam

#index both sorted.bam and sorted_unique.bam files
samtools index ${FILENAME}HISAT2_virome_alignment_Cufflinks_sorted_unique.bam
samtools index ${FILENAME}HISAT2_virome_alignment_Cufflinks_sorted.bam

#count reads in both sorted.bam and sorted_unique.bam files
samtools idxstats ${FILENAME}HISAT2_virome_alignment_Cufflinks_sorted_unique.bam > ${FILENAME}_virus_counts_unique.txt
samtools idxstats ${FILENAME}HISAT2_virome_alignment_Cufflinks_sorted.bam > ${FILENAME}_virus_counts.txt
