#!/bin/bash
#SBATCH --job-name="fastqc"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --output="fastqc_trans2017_2_out.%A-%a"
#SBATCH --error="fastqc_trans2017_2_out.%A-%a"

set -e 
echo "START $(date)"
module load FastQC/0.11.5-Java-1.8.0_92

project_home="/data3/marine_diseases_lab/tejashree/Bio_project_SRA"
sra="$project_home/sra"
qc="$project_home/qc"
#fastqc="$qc/fastqc"
trans2017="$project_home/qc/trans2017_2"
fastqc="$trans2017/fastqc"
source $project_home/github/oyster/lib/slack.sh

cd $fastqc

#files=( _ SRR5357619 SRR5357620 SRR5357623 SRR5357624 SRR5357625 SRR5357626 )
#files=( _ SRR5357623 SRR5357624 SRR5357625 SRR5357626 )
#files=(
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

files=(
    _
    C-K-0_L007
    C-K-0_L008
    C-M-0_L007
    C-M-0_L008
    C-R1_L007
    C-R1_L008
    C-R2_L007
    C-R2_L008
    C-R3_L008
    C-V-0_L007
    C-V-0_L008
    RE-K-6_L007
    RE-K-6_L008
    RE-M-6_L007
    RE-M-6_L008
    RE-R1_L007
    RE-R1_L008
    RE-R2_L007
    RE-R2_L008
    RE-R3_L007
    RE-R3_L008
    RE-V-6_L007
    RE-V-6_L008
    RI-K-24_L007
    RI-K-24_L008
    RI-K-6_L007
    RI-K-6_L008
    RI-M-24_L007
    RI-M-24_L008
    RI-M-6_L007
    RI-M-6_L008
    RIplusRE-R1_L007
    RIplusRE-R1_L008
    RIplusRE-R2_L007
    RIplusRE-R2_L008
    RIplusRE-R3_L007
    RIplusRE-R3_L008
    RI-V-24_L007
    RI-V-24_L008
    RI-V-6_L007
    RI-V-6_L008
    S4-K-24_L007
    S4-K-24_L008
    S4-K-6_L007
    S4-K-6_L008
    S4-M-24_L007
    S4-M-24_L008
    S4-M-6_L007
    S4-M-6_L008
    S4plusRE-R1_L007
    S4plusRE-R1_L008
    S4plusRE-R2_L007
    S4plusRE-R2_L008
    S4plusRE-R3_L007
    S4plusRE-R3_L008
    S4-V-24_L007
    S4-V-24_L008
    S4-V-6_L007
    S4-V-6_L008
    )
# the file prefix for this job is indexing files array by the value of PBS_ARRAYID

me="fastqc::${SLURM_ARRAY_TASK_ID}"

#array=($(ls -1 ${qc}/${files[$PBS_ARRAYID]}_*.fq))
array=($(ls -1 ${trans2017}/${files[$SLURM_ARRAY_TASK_ID]}_*.fq))
echo "[Job $SLURM_ARRAY_TASK_ID] Starting fastqc on ${array[@]}"

for i in ${array[@]}; do  # @ symbol tells it to go through each item in the array  
    post_slack_message cluster-jobs "$me starting Job [$SLURM_ARRAY_TASK_ID] running fastqc ${i}"
    fastqc ${i}
    post_slack_message cluster-jobs "$me Job [$SLURM_ARRAY_TASK_ID] done"
done

post_slack_message cluster-jobs "$me Job [$SLURM_ARRAY_TASK_ID] done"
