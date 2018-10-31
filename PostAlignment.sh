#!/bin/bash
#SBATCH --job-name=PostAlignment   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH -c 1                  # Run a single task
#SBATCH --mem-per-cpu=100gb           # Memory per node
#SBATCH --time=6-23:59:59             # Time limit hrs:min:sec
#SBATCH --output=PostAlignment_%A-%a.out    # Standard output and error log
#SBATCH --array=0-106                 # Array range

#load modules
module purge
module load centos6/0.0.1-fasrc01
module load slurm-drmaa/1.2.0-fasrc01
module load samtools/1.5-fasrc02
module load java/1.8.0_45-fasrc01
module load bedtools2/2.26.0-fasrc01

#list all .sam virome alignment files in directory
FILE=($(ls *HISAT2_virome_alignment_first_pass_Cufflinks.sam  | rev | cut -c 49- | rev | uniq)) 
 
# grab out filename from the array exported from our 'parent' shell
FILENAME=${FILE[$SLURM_ARRAY_TASK_ID]}
echo ${FILENAME}

#convert .sam to .bam
samtools view -Su ${FILENAME}HISAT2_virome_alignment_first_pass_Cufflinks.sam > ${FILENAME}_HISAT2_virome_alignment_first_pass.bam

#sort .bam
samtools sort ${FILENAME}_HISAT2_virome_alignment_first_pass.bam > ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted.bam

#MarkDuplicates with Picard

java -jar /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/GTEx_V7_JQ_Lab/Data/Tools/github.com/broadinstitute/picard/releases/download/2.18.14/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${FILENAME}_HISAT2_virome_alignment_first_pass_sorted.bam O=${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup.bam M=${FILENAME}_marked_dup_metrics.txt 

#capture unique reads
samtools view -bq 50 ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup.bam > ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique.bam

#remove tandem repeats
bedtools intersect -abam ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique.bam -b Shortened_Virome_2.7.7.80.10.50.2000.bed -v > ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique_noTanRep.bam

#index
samtools index ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique_noTanRep.bam
samtools index ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique.bam
samtools index ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup.bam
samtools index ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted.bam

#count reads
samtools idxstats ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique_noTanRep.bam > ${FILENAME}virome_counts_MrkDup_unique_noTanRep.txt
samtools idxstats ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique.bam > ${FILENAME}virome_counts_MrkDup_unique.txt
samtools idxstats ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup.bam > ${FILENAME}virome_counts_MrkDup.txt
samtools idxstats ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted.bam > ${FILENAME}virome_counts.txt
