#!/bin/bash
#SBATCH --job-name=VirusFirstPass   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH -c 1                  # Run a single task
#SBATCH --mem-per-cpu=150gb           # Memory per node
#SBATCH --time=6-23:59:59             # Time limit hrs:min:sec
#SBATCH --output=virus_first_pass_%A-%a.out    # Standard output and error log
#SBATCH --array=0-106                 # Array range

#load modules
module purge
module load centos6/0.0.1-fasrc01
module load slurm-drmaa/1.2.0-fasrc01
module load hisat2/2.1.0-fasrc01

#list all .sam hg38 alignment files in directory
FILE=($(ls *paired.fq |rev | cut -c 21- | rev | uniq))
 
# grab out filename from the array exported from our 'parent' shell
FILENAME=${FILE[$SLURM_ARRAY_TASK_ID]}
echo ${FILENAME}

# run!

hisat2 --phred33 --known-splicesite-infile /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/froy/virome_files/Shortened_Virome_ss.txt --dta-cufflinks  -x /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/froy/virome_files/genome_Virome_index -1 ${FILENAME}___forward_paired.fq -2 ${FILENAME}___reverse_paired.fq -U ${FILENAME}_forward_unpaired.fq,${FILENAME}_reverse_unpaired.fq -S ${FILENAME}HISAT2_virome_alignment_first_pass_Cufflinks.sam --un ${FILENAME}virome_first_pass_Cufflinks_un-seqs --un-conc ${FILENAME}virome_first_pass_Cufflinks_un-conc-mate
