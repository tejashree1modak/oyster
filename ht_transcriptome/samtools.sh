#!/bin/bash
#SBATCH --job-name="samtools"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="samtools_trans2017_final_out.%A-%a"
#SBATCH --error="samtools_trans2017_final_out.%A-%a"

set -e
echo "START" $(date)

module load SAMtools/1.5-foss-2017a

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
bwa=${project_home}/bwa/trans2017_final
#FILES=(
#    _
#    C_K_0_TACAGC
#    C_M_0_CGGAAT
#    C_V_0_CACGAT
#    RE_K_6_TCCCGA
#    RE_M_6_CTCAGA
#    RE_V_6_CATGGC
#    RI_K_24_TCGGCA
#    RI_K_6_TCATTC
#    RI_M_24_TAATCG
#    RI_M_6_CTATAC
#    RI_V_24_CCAACA
#    RI_V_6_CAGGCG
#    S4_K_24_TCGAAG
#    S4_K_6_TATAAT
#    S4_M_24_GACGAC
#    S4_M_6_CTAGCT
#    S4_V_24_CATTTT
#    S4_V_6_CACTCA
#    )

FILES=( 
    _ 
    C_K_0
    C_M_0
    C_V_0
    C_R1
    C_R2
    C_R3
    RE_K_6
    RE_M_6
    RE_V_6
    RE_R1
    RE_R2
    RE_R3
    RI_K_24
    RI_K_6
    RI_M_24
    RI_M_6
    RI_V_24
    RI_V_6
    RIplusRE_R1
    RIplusRE_R2
    RIplusRE_R3
    S4_K_24
    S4_K_6
    S4_M_24
    S4_M_6
    S4_V_24
    S4_V_6
    S4plusRE_R1
    S4plusRE_R2
    S4plusRE_R3
    )

array1=($(ls -1 ${bwa}/${FILES[$SLURM_ARRAY_TASK_ID]}.sam))
echo "Starting flagstat  on ${array1[@]}"
for i in ${array1[@]}; do
    base=$(basename $i)
    samtools flagstat ${i} > ${bwa}/${base%.sam}.stats #get % mapped
    echo " Done '${bwa}/${base%.sam}.stats'"
done
echo "STOP" $(date)
