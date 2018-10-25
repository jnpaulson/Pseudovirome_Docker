#!/bin/bash
#SBATCH --job-name=Trimming   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH -c 1                  # Run a single task
#SBATCH --mem-per-cpu=12gb           # Memory per node
#SBATCH --time=6-23:59:59             # Time limit hrs:min:sec
#SBATCH --output=Trimming_%A-%a.out    # Standard output and error log
#SBATCH --array=0-106                 # Array range

#load modules
module purge
module load centos6/0.0.1-fasrc01
module load slurm-drmaa/1.2.0-fasrc01
module load java/1.8.0_45-fasrc01
module load fastqc/0.11.5-fasrc01

#mkdir for output
mkdir Raw_Reads
mkdir FASTQC_Output

#list all .fq files in directory and print their names
FILE=($(ls *.fastq | rev | cut -c 14- | rev | uniq))

# grab out filename from the array exported from our 'parent' shell
prefix=${FILE[$SLURM_ARRAY_TASK_ID]}
echo ${prefix}

#check on quality of reads pre-trim
./fastqc ${prefix}_R1_001.fastq ${prefix}_R2_001.fastq

# start trimming
java -jar /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/froy/UPenn/www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ${prefix}_R1_001.fastq ${prefix}_R2_001.fastq ${prefix}_forward_paired.fq ${prefix}_forward_unpaired.fq ${prefix}_reverse_paired.fq ${prefix}_reverse_unpaired.fq ILLUMINACLIP:/n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/froy/UPenn/www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36/Merged-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:75

#check on quality of reads post-trim
./fastqc ${prefix}_forward_paired.fq ${prefix}_forward_unpaired.fq ${prefix}_reverse_paired.fq ${prefix}_reverse_unpaired.fq

#mv raw reads and FASTQ output into their respective folders
mv *fastq Raw_Reads
mv *html *zip FASTQC_Output
