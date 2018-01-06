#!/bin/bash
#PBS -l nodes=1
#PBS -l walltime=1000:00:00
#PBS -j oe
#PBS -o preseq_out

set -e
echo "START" $(date)

module load GSL/2.1-foss-2016b

#preseq  code
    #              3. qsub -N name you want to give -t pick a number from 1-4 since our input to FILES array
    #                 has 6 files. eg: 1. RUN_STEP=preseq
    #                                     qsub -N preseq -t 1-4 preseq.sh
    #To check status: qstat -n -1 -t -u tejashree
    # -n: shows node job is running on,-1: single line output, -t: show array details, -u user

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
cvir_genome=${project_home}/cvir_genome
bac_genomes=${project_home}/bac_genomes
hisat2=${project_home}/hisat2
trans2016=${hisat2}/cvir/test
trans2017=${hisat2}/trim5/trans2017     # for trans 2016, the output files are in ${hisat2}/cvir and ${hisat2}/cvir/test
#qc=${project_home}/qc
qc="$project_home/qc/trans2017"
preseq=${project_home}/preseq/preseq

source $project_home/github/oyster/lib/slack.sh

FILES=( _ SRR5357619 SRR5357620 SRR5357623 SRR5357624 SRR5357625 SRR5357626 )
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

RUN_STEP=preseq

me="preseq::${PBS_JOBNAME}::$RUN_STEP on $HOSTNAME"

if [ -z ${RUN_STEP} ] ; then
    # the step to run needs to be specified
    echo "($PBS_ARRAYID) run step not specified.."
    exit 1
fi

if [ "${#FILES[@]}" == "0" ] ; then
    echo "($PBS_ARRAYID) ERROR: No files have been specified, exiting.. : ${#FILES[@]} ${FILES[@]}"
    exit 1
fi

# tell system to call post_slack_message with last 50 lines of output
trap "tail -n 50 preseq_out-${PBS_ARRAYID} | post_slack_message cluster-jobs - '$me::log'" ERR
#-------------------------------------------------------------------------------------------
#Run preseq code on all files in FILES variable
if [ "${RUN_STEP}" == "preseq" ] ; then
    #array1=($(ls -1 ${trans2017}/cvir/${FILES[$PBS_ARRAYID]}.bam))
    array1=($(ls -1 ${trans2016}/${FILES[$PBS_ARRAYID]}.bam))
    if [ "${#array1[@]}" == "0" ]; then
        echo "($PBS_ARRAYID) ERROR: No input files matching '${qc}/${FILES[$PBS_ARRAYID]}_1.fq'"
        post_slack_message cluster-jobs "ERROR: No input files matching '${trans2017}/cvir/${FILES[$PBS_ARRAYID]}.bam'" "$me"
        exit 1
    fi

    #outfile=${preseq}/trans2017/${FILES[$PBS_ARRAYID]}.txt
    outfile=${preseq}/trans2016/${FILES[$PBS_ARRAYID]}.txt
    echo "Starting ${preseq}/preseq c_curve -B -v -o ${outfile} ${array1[@]}" "$me"
    post_slack_message cluster-jobs "Starting ${preseq}/preseq c_curve -B -v -o ${outfile} ${array1[@]}" "$me"
    ${preseq}/preseq c_curve -B -v -o ${outfile} ${array1[@]}

    echo "done"
    post_slack_message cluster-jobs "Done preseq c_curve -B -v -o ${outfile} ${array1[@]}" "$me"
fi

