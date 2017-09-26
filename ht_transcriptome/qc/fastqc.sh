#!/bin/bash
#PBS -N fastqc
#PBS -l nodes=1
#PBS -l walltime=1000:00:00
#PBS -j oe
#PBS -o fastqc_out 
#PBS 1-6
#PBS -t 1-4 

set -e 
echo "START $(date)"
module load FastQC/0.11.5-Java-1.8.0_92

project_home="/data3/marine_diseases_lab/tejashree/Bio_project_SRA"
sra="$project_home/sra"
qc="$project_home/qc"
fastqc="$qc/fastqc"

cd $fastqc

#files=( _ SRR5357619 SRR5357620 SRR5357623 SRR5357624 SRR5357625 SRR5357626 )
files=( _ SRR5357623 SRR5357624 SRR5357625 SRR5357626 )

# the file prefix for this job is indexing files array by the value of PBS_ARRAYID

array=($(ls -1 ${qc}/${files[$PBS_ARRAYID]}_*.fq))
echo "[Job $PBS_ARRAYID] Starting fastqc on ${array[@]}"

for i in ${array[@]}; do  # @ symbol tells it to go through each item in the array  
    echo "[Job $PBS_ARRAYID] running fastqc ${i}"
    fastqc ${i}
    echo "[Job $PBS_ARRAYID] fastqc done"
done

echo "[Job $PBS_ARRAYID] Done fastqc on ${array[@]}"


