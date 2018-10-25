#!/bin/bash
#SBATCH --job-name=Penn_Array   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH -c 1                  # Run a single task
#SBATCH --mem-per-cpu=150gb           # Memory per node
#SBATCH --time=6-23:59:59             # Time limit hrs:min:sec
#SBATCH --output=hg38_array_%A-%a.out    # Standard output and error log
#SBATCH --array=0-106                 # Array range

#load modules
module purge
module load centos6/0.0.1-fasrc01
module load slurm-drmaa/1.2.0-fasrc01
module load hisat2/2.1.0-fasrc01

#list all .fq files in directory and print their names
FILE=($(ls *paired.fq |rev | cut -c 21- | rev | uniq))

# grab out filename from the array exported from our 'parent' shell
FILENAME=${FILE[$SLURM_ARRAY_TASK_ID]}
echo ${FILENAME}

# run!

hisat2 --phred33 --known-splicesite-infile /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/froy/hg38_with_HML2_genome_files/hg38_HML2_genes_051918_ss.txt --dta-cufflinks  -x /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/froy/hg38_with_HML2_genome_files/genome_hg38_HML2_index -1 ${FILENAME}___forward_paired.fq -2 ${FILENAME}___reverse_paired.fq -U ${FILENAME}_forward_unpaired.fq,${FILENAME}_reverse_unpaired.fq -S ${FILENAME}_HISAT2_hg38_HML2_alignment.sam --un ${FILENAME}_HISAT2_hg38_HML2_alignment_un-seqs --un-conc ${FILENAME}_HISAT2_hg38_HML2_alignment_un-conc-mate;
