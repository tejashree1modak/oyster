#!/bin/bash
#PBS -N fastqc
#PBS -l nodes=1
#PBS -l walltime=1000:00:00
#PBS -j oe
#PBS -o fastqc_out 
#PBS 1-18
#PBS -t 1-18 

set -e 
echo "START $(date)"
module load FastQC/0.11.5-Java-1.8.0_92

project_home="/data3/marine_diseases_lab/tejashree/Bio_project_SRA"
sra="$project_home/sra"
qc="$project_home/qc"
#fastqc="$qc/fastqc"
trans2017="$project_home/qc/trans2017"
fastqc="$trans2017/fastqc"

cd $fastqc

#files=( _ SRR5357619 SRR5357620 SRR5357623 SRR5357624 SRR5357625 SRR5357626 )
#files=( _ SRR5357623 SRR5357624 SRR5357625 SRR5357626 )
files=(
    _
    C_K_0_TACAGC
    C_M_0_CGGAAT
    C_V_0_CACGAT
    RE_K_6_TCCCGA
    RE_M_6_CTCAGA
    RE_V_6_CATGGC
    RI_K_24_TCGGCA
    RI_K_6_TCATTC
    RI_M_24_TAATCG
    RI_M_6_CTATAC
    RI_V_24_CCAACA
    RI_V_6_CAGGCG
    S4_K_24_TCGAAG
    S4_K_6_TATAAT
    S4_M_24_GACGAC
    S4_M_6_CTAGCT
    S4_V_24_CATTTT
    S4_V_6_CACTCA
    )
# the file prefix for this job is indexing files array by the value of PBS_ARRAYID

#array=($(ls -1 ${qc}/${files[$PBS_ARRAYID]}_*.fq))
array=($(ls -1 ${trans2017}/${files[$PBS_ARRAYID]}_*.fq))
echo "[Job $PBS_ARRAYID] Starting fastqc on ${array[@]}"

for i in ${array[@]}; do  # @ symbol tells it to go through each item in the array  
    echo "[Job $PBS_ARRAYID] running fastqc ${i}"
    fastqc ${i}
    echo "[Job $PBS_ARRAYID] fastqc done"
done

echo "[Job $PBS_ARRAYID] Done fastqc on ${array[@]}"


