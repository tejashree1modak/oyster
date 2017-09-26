#!/bin/bash
#PBS -l nodes=1
#PBS -l walltime=1000:00:00
#PBS -j oe
#PBS -o hisat_out

set -e
echo "START" $(date)

module load HISAT2/2.0.4-foss-2016b   
module load SAMtools/1.3.1-foss-2016b

#HISAT2 code
#Indexing a reference genome and no annotation file
	#create new directory for the HISAT index  and needed files called pipeline_files, and put the genome inside it
	# copy all reads files into this directory as well to ensure easy access by commands

    #To run script:1. Set RUN_STEP= either hisat or samtools then save and quit file. 
    #              2. qsub -N name you want to give -t pick a number from 1-4 since our input to FILES array
    #                 has 4 files. eg: 1. RUN_STEP=hisat
    #                                     qsub -N 23.fq -t 1 hisat.sh
    #                              eg: 2 RUN_STEP=samtools
    #                                     qsub -N 24.sam -t 2 hisat.sh
    #To check status: qstat -n -1 -t -u tejashree
    # -n: shows node job is running on,-1: single line output, -t: show array details, -u user

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
cvir_genome=${project_home}/cvir_genome
hisat2=${project_home}/hisat2
qc=${project_home}/qc

cd $hisat2
mkdir -p cvir

if [ -z ${PBS_ARRAYID} ] ; then
    # this has to run in ARRAY mode
    echo "Run with qsub -t <array args>, exiting.."
    exit 1
fi

source $project_home/github/oyster/lib/slack.sh

FILES=( _ SRR5357619 SRR5357620 SRR5357623 SRR5357624 SRR5357625 SRR5357626 )
#RUN_STEP=hisat
RUN_STEP=samtools
me="hisat::${PBS_JOBNAME}::$RUN_STEP"

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
trap "tail -n 50 hisat_out-${PBS_ARRAYID} | post_slack_message cluster-jobs - '$me::log'" ERR

# ------------------------------------------------------------------------------------------
#hisat2 align paired end reads to reference genome
if [ "${RUN_STEP}" == "hisat" ] ; then
    array1=($(ls -1 ${qc}/${FILES[$PBS_ARRAYID]}_1.fq))
    if [ "${#array1[@]}" == "0" ]; then
        echo "($PBS_ARRAYID) ERROR: No input files matching '${qc}/${FILES[$PBS_ARRAYID]}_1.fq'"
        post_slack_message cluster-jobs "ERROR: No input files matching '${qc}/${FILES[$PBS_ARRAYID]}_1.fq'" "$me"
        exit 1
    fi

    echo "Starting hisat2 alignment on ${array1[@]}"
    post_slack_message cluster-jobs "Starting hisat2 alignment on ${array1[@]}" "$me"
    for left_file in ${array1[@]}; do
        right_file=$(echo ${left_file}|sed s/_1/_2/)
        file=$(basename ${left_file%_1.fq}) 

        echo starting hisat2 --dta -x ${cvir_genome}/cvir -1 ${left_file} -2 ${right_file} -S ${hisat2}/cvir/${file}.sam
        post_slack_message cluster-jobs "starting hisat2 --dta -x ${cvir_genome}/cvir -1 ${left_file} -2 ${right_file} -S ${hisat2}/cvir/${file}.sam" "${me}.${file}"
        if [ "$debug" ] ; then
            touch ${hisat2}/cvir/${file}.sam
        else
            hisat2 --dta -x ${cvir_genome}/cvir -1 ${left_file} -2 ${right_file} -S ${hisat2}/cvir/${file}.sam
        fi
        echo done
        post_slack_message cluster-jobs "done" "${me}.${file}"
    done
    echo "Done hisat2 alignment"
    post_slack_message cluster-jobs "Done hisat2 alignment" "$me"
fi

# ------------------------------------------------------------------------------------------
#SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie
if [ "${RUN_STEP}" == "samtools" ] ; then
    array3=($(ls -1 ${hisat2}/cvir/${FILES[$PBS_ARRAYID]}.sam))
    if [ "${#array3[@]}" == "0" ]; then
        echo "($PBS_ARRAYID) ERROR: No input files matching '${qc}/${FILES[$PBS_ARRAYID]}.sam'"
        post_slack_message cluster-jobs "ERROR: No input files matching '${qc}/${FILES[$PBS_ARRAYID]}.sam'" "$me"
        exit 1
    fi

    echo "Starting samtools convert on ${array3[@]}"
    post_slack_message cluster-jobs "Starting samtools convert on ${array3[@]}" "$me"
    for i in ${array3[@]}; do
        echo starting samtools sort -m 10G -o ${i%.sam}.bam ${i}
        post_slack_message cluster-jobs "starting samtools sort -m 10G -o ${i%.sam}.bam ${i}" "$me"
        if [ "$debug" ] ; then
            touch ${i%.sam}.bam
        else
            samtools sort -o ${i%.sam}.bam ${i}
        fi
        echo done
        post_slack_message cluster-jobs "done" "$me"
    done
    echo "Done samtools convert"
    post_slack_message cluster-jobs "Done samtools conver" "$me"
fi

echo "STOP" $(date)
