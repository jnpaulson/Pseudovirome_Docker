#!/bin/bash
#SBATCH --job-name=PicardDedup   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH -c 1                  # Run a single task
#SBATCH --mem-per-cpu=10gb           # Memory per node
#SBATCH --time=6-23:59:59             # Time limit hrs:min:sec
#SBATCH --output=Picard_Dedup_%A-%a.out    # Standard output and error log
#SBATCH --array=0-106                 # Array range

#load modules
module purge
module load centos6/0.0.1-fasrc01
module load slurm-drmaa/1.2.0-fasrc01
module load java/1.8.0_45-fasrc01

#list all .fq files in directory and print their names
FILE=($(ls *HISAT2_virome_alignment_Cufflinks_sorted.bam |rev | cut -c 45- | rev | uniq))

# grab out filename from the array exported from our 'parent' shell
FILENAME=${FILE[$SLURM_ARRAY_TASK_ID]}
echo ${FILENAME}

# run!

java -jar /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/GTEx_V7_JQ_Lab/Data/Tools/github.com/broadinstitute/picard/releases/download/2.18.14/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted.bam O=${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup.bam M=${FILENAME}_marked_dup_metrics.txt 
