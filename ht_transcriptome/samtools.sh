#!/bin/bash
#PBS -l nodes=1
#PBS -l walltime=1000:00:00
#PBS -j oe
#PBS -o samtools_out

set -e
echo "START" $(date)

module load SAMtools/1.5-foss-2017a

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
bwa=${project_home}/bwa/trans2017
FILES=(
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

array1=($(ls -1 ${bwa}/${FILES[$PBS_ARRAYID]}.sam))
echo "Starting flagstat  on ${array1[@]}"
for i in ${array1[@]}; do
    base=$(basename $i)
    samtools flagstat ${i} > ${bwa}/${base%.sam}.stats #get % mapped
    echo " Done '${bwa}/${base%.sam}.stats'"
done
echo "STOP" $(date)
