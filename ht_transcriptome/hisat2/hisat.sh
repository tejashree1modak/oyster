#!/bin/bash
#SBATCH --job-name="hisat"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="hisat_trans2017_final_out.%A-%a"
#SBATCH --error="hisat_trans2017_final_out.%A-%a"

set -e
echo "START" $(date)

###
# Run the script like this - 
#
# sbatch -a 5-11,13,15,16,18-22%10 --ntasks-per-node=1 -c 20 /data3/marine_diseases_lab/tejashree/Bio_project_SRA/github/oyster/ht_transcriptome/hisat2/hisat.sh
# 
# ntasks-per-node is critical

module load HISAT2/2.0.4-foss-2016b   
module load SAMtools/1.3.1-foss-2016b

#HISAT2 code
#Indexing a reference genome and no annotation file
    #To run script:1. Set RUN_STEP= either hisat or samtools then save and quit file.
    #              2. If making genome index just set RUN_STEP=genome_index and qsub -N genomeindex hisat.sh
    #              3. qsub -N name you want to give -t pick a number from 1-4 since our input to FILES array
    #                 has 4 files. eg: 1. RUN_STEP=hisat

    #                                     qsub -N 24.sam -t 2 hisat.sh
    #To check status: qstat -n -1 -t -u tejashree
    # -n: shows node job is running on,-1: single line output, -t: show array details, -u user

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
cvir_genome=${project_home}/cvir_genome
bac_genomes=${project_home}/bac_genomes
hisat2=${project_home}/hisat2
trans2017=${hisat2}/trim5/trans2017_final/     # for trans 2016, the output files are in ${hisat2}/cvir and ${hisat2}/cvir/test
#qc=${project_home}/qc
qc="$project_home/qc/trans2017_final"

#cd $hisat2
#mkdir -p cvir
source $project_home/github/oyster/lib/slack.sh

#FILES=( _ SRR5357619 SRR5357620 SRR5357623 SRR5357624 SRR5357625 SRR5357626 )
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

#RUN_STEP=hisat
#RUN_STEP=samtools
#RUN_STEP=samtools_stats
#RUN_STEP=samtools_index

me="hisat::${SBATCH_JOB_NAME}-${SLURM_ARRAY_TASK_ID}::$RUN_STEP on $HOSTNAME"

if [ -z ${RUN_STEP} ] ; then
    # the step to run needs to be specified
    echo "($SLURM_ARRAY_TASK_ID) run step not specified.."
    exit 1
fi

if [ "${#FILES[@]}" == "0" ] ; then
    echo "($SLURM_ARRAY_TASK_ID) ERROR: No files have been specified, exiting.. : ${#FILES[@]} ${FILES[@]}"
    exit 1
fi

# tell system to call post_slack_message with last 50 lines of output
trap "tail -n 50 hisat_trans2017_final_out.*-${SLURM_ARRAY_TASK_ID} | post_slack_message cluster-jobs - '$me::log'" ERR
#-------------------------------------------------------------------------------------------
#Indexing a reference genome and no annotation file
	#create new directory for the HISAT index  and needed files called pipeline_files, and put the genome inside it
	# copy all reads files into this directory as well to ensure easy access by commandsi
#if [ "${RUN_STEP}" == "genome_index" ] ; then
#    cd ${cvir_genome}
#    hisat2-build -f ${cvir_genome}/cvir_edited.fa  cvir_edited_index
#    echo "genome index built"
#fi
if [ "${RUN_STEP}" == "genome_index" ] ; then
    cd ${bac_genomes}
    hisat2-build -f ${bac_genomes}/ALLgenomes_bac.fna  allgenomes_bac_index
    echo "genome index built"
fi

    # -f indicates that the reference input files are FASTA files

    #Stay in the directory created in the previous step

# ------------------------------------------------------------------------------------------
#hisat2 align paired end reads to reference genome
if [ "${RUN_STEP}" == "hisat" ] ; then
    array1=($(ls -1 ${qc}/${FILES[$SLURM_ARRAY_TASK_ID]}_1.fq))
    if [ "${#array1[@]}" == "0" ]; then
        echo "($SLURM_ARRAY_TASK_ID) ERROR: No input files matching '${qc}/${FILES[$SLURM_ARRAY_TASK_ID]}_1.fq'"
        post_slack_message cluster-jobs "ERROR: No input files matching '${qc}/${FILES[$SLURM_ARRAY_TASK_ID]}_1.fq'" "$me"
        exit 1
    fi

    echo "Starting hisat2 alignment on ${array1[@]}"
    post_slack_message cluster-jobs "Starting hisat2 alignment on ${array1[@]}" "$me"
    for left_file in ${array1[@]}; do
        right_file=$(echo ${left_file}|sed s/_1/_2/)
        file=$(basename ${left_file%_1.fq}) 

        #echo starting hisat2 --dta -x ${cvir_genome}/cvir_edited_index -1 ${left_file} -2 ${right_file} -S ${trans2017}/cvir/${file}.sam
        echo starting hisat2 --dta -x ${bac_genomes}/allgenomes_bac_index -1 ${left_file} -2 ${right_file} -S ${trans2017}/bac/${file}.sam
        #post_slack_message cluster-jobs "starting hisat2 --dta -x ${cvir_genome}/cvir_edited_index -1 ${left_file} -2 ${right_file} -S ${trans2017}/cvir/${file}.sam" "${me}.${file}"
        post_slack_message cluster-jobs "starting hisat2 --dta -x ${bac_genomes}/allgenomes_bac_index -1 ${left_file} -2 ${right_file} -S ${trans2017}/bac/${file}.sam" "${me}.${file}"
        if [ "$debug" ] ; then
            touch ${trans2017}/bac/${file}.sam
        else
            #hisat2 --dta -x ${cvir_genome}/cvir_edited_index -1 ${left_file} -2 ${right_file} -S ${trans2017}/cvir/${file}.sam
            hisat2 --dta -x ${bac_genomes}/allgenomes_bac_index -1 ${left_file} -2 ${right_file} -S ${trans2017}/bac/${file}.sam
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
    array3=($(ls -1 ${trans2017}/cvir/${FILES[$SLURM_ARRAY_TASK_ID]}.sam))
    if [ "${#array3[@]}" == "0" ]; then
        echo "($SLURM_ARRAY_TASK_ID) ERROR: No input files matching '${trans2017}/${FILES[$SLURM_ARRAY_TASK_ID]}.sam'"
        post_slack_message cluster-jobs "ERROR: No input files matching '${trans2017}/${FILES[$SLURM_ARRAY_TASK_ID]}.sam'" "$me"
        exit 1
    fi

    echo "Starting samtools convert on ${array3[@]}"
    post_slack_message cluster-jobs "Starting samtools convert on ${array3[@]}" "$me"
    for i in ${array3[@]}; do
        echo starting samtools sort -m 50G -o ${i%.sam}.bam ${i}
        post_slack_message cluster-jobs "starting samtools sort -m 50G -o ${i%.sam}.bam ${i}" "$me"
        if [ "$debug" ] ; then
            touch ${i%.sam}.bam
        else
            samtools sort -m 50G -o ${i%.sam}.bam ${i}
        fi
        echo done
        post_slack_message cluster-jobs "done" "$me"
    done
    echo "Done samtools convert"
    post_slack_message cluster-jobs "Done samtools conver" "$me"
fi
#-----------------------------------------------------------------------------
#Get bam file statistics for percentage aligned with flagstat
# to get more detailed statistics use $ samtools stats ${i}

if [ "${RUN_STEP}" == "samtools_stats" ] ; then
    array3=($(ls -1 ${trans2017}/${FILES[$SLURM_ARRAY_TASK_ID]}.bam))
    for i in ${array3[@]}; do
        base=$(basename $i)
        samtools flagstat ${i} > ${trans2017}/${base%.bam}.stats #get % mapped
    #to extract more detailed summary numbers
        samtools stats ${i}| grep ^SN | cut -f 2- > ${trans2017}/${base%.bam}.fullstat
  #| grep ^SN | cut -f 2- > ${trans2017}/${i}.bam.fullstat
        echo "STATS DONE" $(date)
    done
fi
#------------------------------------------------------------------
#samtools index bam files 
if [ "${RUN_STEP}" == "samtools_index" ] ; then
    array3=($(ls -1 ${trans2017}/${FILES[$SLURM_ARRAY_TASK_ID]}.bam))
    for i in ${array3[@]}; do
        base=$(basename $i)
        samtools index -b ${i} ${trans2017}/${base%.bam}.index 
        echo "index done" $(date)
    done
fi
echo "STOP" $(date)
