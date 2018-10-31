#!/bin/bash
#SBATCH --job-name=Removing_hg38_Reads   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH -c 1                  # Run a single task
#SBATCH --mem-per-cpu=100gb           # Memory per node
#SBATCH --time=6-23:59:59             # Time limit hrs:min:sec
#SBATCH --output=Removing_hg38_Reads_%A-%a.out    # Standard output and error log
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

#get read names from sam file
cut -f1 ${FILENAME}_HISAT2_hg38_HML2_alignment_second_pass.sam > ${FILENAME}HISAT2_hg38_HML2_alignment_second_pass_readnames.txt

#remove reads in text file from filtered virome bam file
java -jar /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/GTEx_V7_JQ_Lab/Data/Tools/github.com/broadinstitute/picard/releases/download/2.18.14/picard.jar FilterSamReads I=${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique_noTanRep.bam O=${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique_noTanRep_hg38readsremoved.bam READ_LIST_FILE=${FILENAME}HISAT2_hg38_HML2_alignment_second_pass_readnames.txt FILTER=excludeReadList 

#index
samtools index ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique_noTanRep_hg38readsremoved.bam

#count reads
samtools idxstats ${FILENAME}_HISAT2_virome_alignment_first_pass_sorted_MrkDup_unique_noTanRep_hg38readsremoved.bam > ${FILENAME}virome_counts_MrkDup_unique_noTanRep_hg38readsremoved.txt
