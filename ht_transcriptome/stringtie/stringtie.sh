#!/bin/bash
#SBATCH --job-name="stringtie"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="stringtie_trans2017_final_out.%A-%a"
#SBATCH --error="stringtie_trans2017_final_out.%A-%a"

set -e
echo "START" $(date)

#This script takes bam files from HISAT (processed by SAMtools) and performs StringTie assembly and quantification and converts
# data into a format that is readable as count tables for DESeq2 usage
#To run script:1. Set RUN_STEP= either assembly then save and quit file. 
    #              2. qsub -N name you want to give -t pick a number from 1-4 since our input to FILES array
    #                 has 4 files. eg: 1. RUN_STEP=assembly
    #                                     qsub -N 23.bam -t 1 stringtie.sh
    #                              eg: 2 RUN_STEP=assembly
    #                                     qsub -N 24.bam -t 2 stringtie.sh
    #To check status: qstat -n -1 -t -u tejashree
    # -n: shows node job is running on,-1: single line output, -t: show array details, -u user

module load StringTie/1.3.3b-foss-2016b
module load gffcompare/0.10.1-foss-2016b
echo "modules loaded"

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA
stringtie=${project_home}/stringtie
hisat2_out=${project_home}/hisat2/cvir/test/
cvir_genome=${project_home}/cvir_genome
hisat2=${project_home}/hisat2
hisat_cvir=${hisat2}/cvir
bam_files=${project_home}/hisat2/trim5/trans2017_final/cvir
cd ${stringtie}/test/trans2017_final

source $project_home/github/oyster/lib/slack.sh


#FILES=( _ SRR5357619 SRR5357620 SRR5357623 SRR5357624 SRR5357625 SRR5357626 )
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

#RUN_STEP=assembly
#RUN_STEP=merge
#RUN_STEP=compare
#RUN_STEP=reestimate
#RUN_STEP=deseq

me="stringtie::${SBATCH_JOB_NAME}-${SLURM_ARRAY_TASK_ID}::$RUN_STEP on $HOSTNAME"

if [ -z ${RUN_STEP} ] ; then
    # the step to run needs to be specified
    echo "($SLURM_ARRAY_TASK_ID) run step not specified.."
    exit 1
fi

# tell system to call post_slack_message
trap "tail -n 50 stringtie_out-${SLURM_ARRAY_TASK_ID} | post_slack_message cluster-jobs - '$me::log'" ERR

## ------------------------------------------------------------------------------------------
# StringTie to assemble transcripts for each sample with the GFF3 annotation file
if [ "${RUN_STEP}" == "assembly" ] ; then
    #array1=($(ls -1 ${hisat_cvir}/test/${FILES[$SLURM_ARRAY_TASK_ID]}.bam))
    array1=($(ls -1 ${bam_files}/${FILES[$SLURM_ARRAY_TASK_ID]}.bam))
    if [ "${#array1[@]}" == "0" ]; then
        echo "($SLURM_ARRAY_TASK_ID) ERROR: No input files matching '${bam_files}/${FILES[$SLURM_ARRAY_TASK_ID]}.bam'"
        post_slack_message cluster-jobs "ERROR: No input files matching '${bam_files}/${FILES[$SLURM_ARRAY_TASK_ID]}.bam'" "$me"
        exit 1
    fi

    echo "Starting assembly"
    post_slack_message cluster-jobs "Starting assembly on ${array1[@]}" "$me"
    for i in ${array1[@]}; do
        base=$(basename $i)
        echo stringtie -G ${cvir_genome}/ref_C_virginica-3.0_top_level.gff3 \
                  -o $stringtie/test/trans2017_final/${base%.bam}.gtf -p 10 -l ${base%.bam} ${i}
        post_slack_message cluster-jobs "stringtie -G ${cvir_genome}/ref_C_virginica-3.0_top_level.gff3 -o $stringtie/test/trans2017_final/${base%.bam}.gtf -p 10 -l ${base%.bam} ${i}" "$me"
        if [ "$debug" ]; then
            touch ${i%.bam}.gtf 
        else
            #stringtie -G ${cvir_genome}/ref_C_virginica-3.0_top_level.gff3 \
            #  -o $stringtie/test/${base%.bam}.gtf -l ${base%.bam} -p 10 ${i}
            stringtie -G ${cvir_genome}/ref_C_virginica-3.0_top_level.gff3 \
              -o $stringtie/test/trans2017_final/${base%.bam}.gtf -l ${base%.bam} -p 10 ${i}
        fi
        echo done
        post_slack_message cluster-jobs "done" "$me.$base"
    done 
    post_slack_message cluster-jobs "Done assembly on ${array1[@]}" "$me"
fi

	# command structure: $ stringtie <options> -G <reference.gtf or .gff> -o outputname.gtf -l prefix_for_transcripts input_filename.bam
	# -o specifies the output name
	# -G specifies you are aligning with an option GFF or GTF file as well to perform novel transcript discovery 
	# -l Sets <label> as the prefix for the name of the output transcripts. Default: STRG
	# don't use -e here if you want it to assemble any novel transcripts
#----------------------------------------------------------------------------------------------------	
#StringTie Merge, will merge all GFF files and assemble transcripts into a non-redundant set of transcripts, after which re-run StringTie with -e
	
	#create mergelist.txt in nano, names of all the GTF files created in the last step with each on its own line
	#ls *.gtf > mergelist.txt

	#check to sure one file per line
	#cat mergelist.txt

	#Run StringTie merge, merge transcripts from all samples (across all experiments, not just for a single experiment)
if [ "${RUN_STEP}" == "merge" ] ; then
    echo "Starting StringTie Merge"
    post_slack_message cluster-jobs "Starting merge" "$me"
 	stringtie --merge -G ${cvir_genome}/ref_C_virginica-3.0_top_level.gff3 -o $stringtie/test/trans2017_final/stringtie_merged.gtf $stringtie/test/trans2017_final/mergelist.txt
    echo "done"
    post_slack_message cluster-jobs "done" "$me.$base"
    #done
fi
#-----------------------------------------------------------------------
#Cuffcompare to 
#gffcompare to compare how transcripts compare to reference annotation
if [ "${RUN_STEP}" == "compare" ] ; then
    echo "Starting StringTie Gffcompare"
    post_slack_message cluster-jobs "Starting gffcompare" "$me"
 	gffcompare -r ${cvir_genome}/ref_C_virginica-3.0_top_level.gff3 -G -o $stringtie/merged $stringtie/test/trans2017_final/stringtie_merged.gtf
    echo "done"
    post_slack_message cluster-jobs "compare done" "$me.$base"
    #done
fi
    # -o specifies prefix to use for output files
	# -r followed by the annotation file to use as a reference
 	# merged.annotation.gtf tells you how well the predicted transcripts track to the reference annotation file
 	# merged.stats file shows the sensitivity and precision statistics and total number for different features (genes, exons, transcripts)
#----------------------------------------------------------------------------
#Re-estimate transcript abundance after merge step

if [ "${RUN_STEP}" == "reestimate" ] ; then
    array1=($(ls -1 ${bam_files}/${FILES[$SLURM_ARRAY_TASK_ID]}.bam))
    if [ "${#array1[@]}" == "0" ]; then
        echo "($SLURM_ARRAY_TASK_ID) ERROR: No input files matching '${bam_files}/${FILES[$SLURM_ARRAY_TASK_ID]}.bam'"
        post_slack_message cluster-jobs "ERROR: No input files matching '${bam_files}/${FILES[$SLURM_ARRAY_TASK_ID]}.bam'" "$me"
        exit 1
    fi
    
    echo "Starting re-estimation"
    post_slack_message cluster-jobs "Starting re-estimation on ${array1[@]}" "$me"
    for i in ${array1[@]}; do
        echo "${array1[@]}"
        base=$(basename $i)
        echo stringtie -e -G $stringtie/test/trans2017_final/stringtie_merged.gtf \
                  -o $stringtie/test/trans2017_final/${i%.bam}.merge.gtf -p 10 ${i}
        post_slack_message cluster-jobs "stringtie -e -G $stringtie/test/trans2017_final/stringtie_merged.gtf -o $stringtie/test/trans2017_final/${base%.bam}.merge.gtf -p 10 ${i}" "$me"
        if [ "$debug" ]; then
            touch ${i%.bam}.merge.gtf 
        else
            stringtie -e -G $stringtie/test/trans2017_final/stringtie_merged.gtf \
              -o $stringtie/test/trans2017_final/${base%.bam}.merge.gtf -p 10 ${i}
        fi
        echo done
        post_slack_message cluster-jobs "done" "$me.$base"
    done 
    post_slack_message cluster-jobs "Done re-estimation on ${array1[@]}" "$me"
fi    
 #erin   
#    for i in ${array1[@]}; do
#		stringtie -e -G $stringtie/stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
#		echo "${i}"
#	done 
	# input here is the original set of alignment files
	# here -G refers to the merged GTF files
	# -e creates more accurate abundance estimations with input transcripts, needed when converting to     #DESeq2 tables
#-----------------------------------------------------------------------------

# Protocol to generate count matrices for genes and transcripts for import into DESeq2 using (prepDE.py) to extract this read count information directly from the files generated by StringTie (run with the -e parameter).
	#generates two CSV files containing the count matrices for genes and transcripts, Given a list of GTFs, which were re-estimated upon merging create sample_list.txt

#Generate count matrices using prepDE.py, prep_DE.py accepts a .txt file listing sample IDs and GTFs paths 
#create sample_list.txt
if [ "${RUN_STEP}" == "deseq" ] ; then
    array2=($(ls $stringtie/test/trans2017_final/*.merge.gtf))
    :> sample_list.txt #to avoid concatenating file content if run twice this command will 
                       #make a new file each time you run the code.
    
    for i in ${array2[@]}; do
        base=$(basename $i)
    	echo "${base%.merge.gtf} ${i}" >> sample_list.txt
    done
    python $stringtie/prepDE.py -i sample_list.txt
fi    	
			
#Steps continued in R, 08_DESeq_RNA_pipeline_just_pathogen_challenge.R


#Helpful reference website http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
