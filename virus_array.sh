#!/bin/bash
#SBATCH --job-name=Virus_PA   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH -c 1                  # Run a single task
#SBATCH --mem-per-cpu=10gb           # Memory per node
#SBATCH --time=6-23:59:59             # Time limit hrs:min:sec
#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH --array=0-89%30                 # Array range

#load modules
module purge
module load centos6/0.0.1-fasrc01
module load slurm-drmaa/1.2.0-fasrc01
module load hisat2/2.1.0-fasrc01

#list all .fq files in directory and print their names
FILE=($(ls *HISAT2_hg38_HML2_alignment_combined.fq |rev | cut -c 39- | rev | uniq))
echo ${FILE}

# grab out filename from the array exported from our 'parent' shell
FILENAME=${FILE[$SLURM_ARRAY_TASK_ID]}
echo ${FILENAME}

# run!

hisat2 --phred33 --known-splicesite-infile /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/froy/virome_files/Shortened_Virome_ss.txt --dta-cufflinks  -x /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/froy/virome_files/genome_Virome_index  -q ${FILENAME}HISAT2_hg38_HML2_alignment_combined.fq -S ${FILENAME}HISAT2_virome_alignment_Cufflinks.sam --un ${FILENAME}virome_Cufflinks_un-seqs --un-conc ${FILENAME}virome_Cufflinks_un-conc-mate