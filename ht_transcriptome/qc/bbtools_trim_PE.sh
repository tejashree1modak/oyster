#!/bin/bash
#PBS-l nodes=1
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -o bbtools_out

# This script processes SRA PE end reads with BBtools to find adaptor sequences
#  with BBmerge, and then uses these for adaptor trimming and quality trimming with bbduk.sh.

set -e 
echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

project_home="/data3/marine_diseases_lab/tejashree/Bio_project_SRA"
sra="$project_home/trans2017/Fastq"
qc="$project_home/qc/trans2017"
source $project_home/github/oyster/lib/slack.sh

#all files are in the home directory and either have ending 
# _1.fq or _2.fq
#going to make two array variables and then iterate through them as an index
#changing all the file endings to .fq to see if bbduk.sh prefers that file name ending

#RUN_STEP=stat
#RUN_STEP=adapter-trim
RUN_STEP=quality-filter

me="bbtools::${PBS_JOBNAME}-$RUN_STEP"

if [ -z ${RUN_STEP} ] ; then
    # the step to run needs to be specified
    echo "($PBS_ARRAYID) run step not specified.."
    exit 1
fi

# tell system to call post_slack_message
trap "tail -n 50 stringtie_out-${PBS_ARRAYID} | post_slack_message cluster-jobs - '$me::log'" ERR

if [ "${RUN_STEP}" == "stat" ]; then
    left=($(ls -1 $sra/*_1.fastq))
    if [ -n "$PBS_ARRAYID" ]; then
        ix=$PBS_ARRAYID
        let ix--
        left=( ${left[$ix]} )
    fi

    echo "Starting 'stat'"
    for left_file in ${left[@]}; do  # @ symbol tells it to go through each item in the array  
        left_file_basename=$(basename $left_file)
        right_file=$(echo ${left_file}|sed s/_1/_2/)
        file=${left_file_basename%_1.fastq} 

        echo "bbduk.sh in1=${left_file} in2=${right_file} k=23 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa stats=${qc}/${file}.stat out=${qc}/${file}.out"
        if [ "$debug" ] ; then
            touch ${qc}/${file}.stat ${qc}/${file}.out 
        else
            post_slack_message cluster-jobs "Starting stat on $left_file  $right_file $me"
            bbduk.sh in1=${left_file} in2=${right_file} k=23 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa stats=${qc}/${file}.stat out=${qc}/${file}.out
        fi
        echo done
    done
    post_slack_message cluster-jobs "Done state $me"
    echo "Done 'stat'"
fi

# stats.txt will then list the names of adapter sequences found, and their frequency
##need to find out if the adaptor sequences are on the left side or the right side
#
#echo "Starting 'Force-trim'"
#for left_file in ${left[@]}; do  # @ symbol tells it to go through each item in the array  
#    left_file_basename=$(basename $left_file)
#    right_file=$(echo ${left_file}|sed s/_1/_2/)
#
#    #Force-trim modulo:clips last base not removed by Illumina (should be 2x150)
#    echo starting bbduk.sh in1=${left_file}  out1=${qc}/${left_file_basename%.fastq}.ftm \
#                 in2=${right_file} out2=${qc}/${left_file_basename%_1.fastq}_2.ftm ftm=5
#    if [ "$debug" ] ; then
#        touch ${qc}/${left_file_basename%.fastq}.ftm ${qc}/${left_file_basename%_1.fastq}_2.ftm 
#    else
#        bbduk.sh in1=${left_file}  out1=${qc}/${left_file_basename%.fastq}.ftm \
#                 in2=${right_file} out2=${qc}/${left_file_basename%_1.fastq}_2.ftm ftm=5
#    fi
#    echo done
#done
#echo "Done 'Force-trim'"
#
## Trimming of adaptors found in the previous command

if [ "${RUN_STEP}" == "adapter-trim" ]; then
    array3=($(ls -1 ${sra}/*_1.fastq))
    if [ -n "$PBS_ARRAYID" ]; then
        ix=$PBS_ARRAYID
        let ix--
        array3=( ${array3[$ix]} )
    fi
    echo "Starting adaptor trim"
    for left_file in ${array3[@]}; do  # @ symbol tells it to go through each item in the array
        right_file=$(echo ${left_file}|sed s/_1/_2/)
        echo starting bbduk.sh in1=${left_file}  out1=${left_file%.fastq}.adp \
                 in2=${right_file} out2=${right_file%.fastq}.adp \
                 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa \
                 ktrim=r k=23 mink=11 hdist=1 tpe tbo
        post_slack_message cluster-jobs "$me: starting bbduk.sh in1=${left_file}  out1=${left_file%.fastq}.adp in2=${right_file} out2=${right_file%.fastq}.adp ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo"
        if [ "$debug" ]; then
            touch ${left_file%.fastq}.adp ${right_file%.fastq}.adp 
        else
            bbduk.sh in1=${left_file}  out1=${left_file%.fastq}.adp \
                 in2=${right_file} out2=${right_file%.fastq}.adp \
                 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa \
                 ktrim=r k=23 mink=11 hdist=1 tpe tbo
        fi
        echo done
    done	
    echo "Done adaptor trim"
fi

#ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
#hdist = hamming distance, hdist =1 allows for 1 mismatch
#flag -tbo specifies to also trim adaptors based on pir overlap detection using BBMerge 
#which does not require known adapter sequences)
#flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)


##quality trimming, of both the left and the right sides to get rid of reads that are less than quality 10
#echo "Starting quality trimming"
#array4=($(ls -1 ${qc}/*_1.adp))	
#for left_file in ${array4[@]}; do
#    right_file=$(echo ${left_file}|sed s/_1/_2/)
#    echo starting bbduk.sh in1=${left_file}  out1=${left_file%.adp}.trim \
#             in2=${right_file} out2=${right_file%.adp}.trim qtrim=rl trimq=20
#
#    if [ "$debug" ]; then
#        touch ${left_file%.adp}.trim ${right_file%.adp}.trim 
#    else
#        bbduk.sh in1=${left_file}  out1=${left_file%.adp}.trim \
#             in2=${right_file} out2=${right_file%.adp}.trim qtrim=rl trimq=20
#    fi
#done
#echo "Done quality trimming"

#
#quality filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10
#echo "Starting quality filtering"
#

# ------------------------------------------------------------------------------------------

if [ "${RUN_STEP}" == "quality-filter" ]; then
    ##array5=($(ls -1 ${qc}/*_1.trim))
    #array5=( ${qc}/SRR5357626_1.trim  )
    array5=($(ls -1 ${qc}/*_1.adp))
    if [ -n "$PBS_ARRAYID" ]; then
        ix=$PBS_ARRAYID
        let ix--
        array5=( ${array5[$ix]} )
    fi
    
    for left_file in ${array5[@]}; do 
        right_file=$(echo ${left_file}|sed s/_1/_2/)
        echo starting bbduk.sh in1=${left_file}  out1=${left_file%.adp}.ftl \
                 in2=${right_file} out2=${right_file%.adp}.ftl maq=10
        post_slack_message cluster-jobs "$me: starting bbduk.sh in1=${left_file}  out1=${left_file%.adp}.ftl in2=${right_file} out2=${right_file%.adp}.ftl maq=10"
         
        if [ "$debug" ]; then
            touch ${left_file%.adp}.ftl ${right_file%.adp}.ftl 
        else
            bbduk.sh in1=${left_file}  out1=${left_file%.adp}.ftl \
                 in2=${right_file} out2=${right_file%.adp}.ftl maq=10
        fi
        echo done
        post_slack_message cluster-jobs "$me: $left_file $right_file done"
    done
    post_slack_message cluster-jobs "$me: done quality filtering"
    echo "Done quality filtering"
fi

#trim first 10 bases 
#echo "Starting trim first 10"
#array6=($(ls -1 ${qc}/*_1.ftl))
#array6=($(ls -1 ${qc}/SRR5357619_1.ftl))
#for left_file in ${array6[@]}; do
#    right_file=$(echo ${left_file}|sed s/_1/_2/)
#    echo starting bbduk.sh in1=${left_file}  out1=${left_file%.ftl}.ftm10 \
#             in2=${right_file} out2=${right_file%.ftl}.ftm10 ftl=10
#    if [ "$debug" ]; then
#        touch ${left_file%.ftl}.ftm10 ${right_file%.ftl}.ftm10 
#    else
#        bbduk.sh in1=${left_file}  out1=${left_file%.ftl}.ftm10 \
#             in2=${right_file} out2=${right_file%.ftl}.ftm10 ftl=10
#    fi
#    echo done
#done 
#echo "Done trim first 10"
##
#echo "Starting remove read less than 25bp"
##remove reads less than 25bp long 
#array7=($(ls -1 ${qc}/*_1.ftm10))
#array7=($(ls -1 ${qc}/SRR5357626_1.ftm10))
#for left_file in ${array7[@]}; do
#    right_file=$(echo ${left_file}|sed s/_1/_2/)
#    echo bbduk.sh in1=${left_file} out1=${left_file%.ftm10}.fq \
#         in2=${right_file} out2=${right_file%.ftm10}.fq minlen=25
#    if [ "$debug" ]; then
#        touch ${left_file%.ftm10}.fq ${right_file%.ftm10}.fq 
#    else
#        bbduk.sh in1=${left_file} out1=${left_file%.ftm10}.fq \
#             in2=${right_file} out2=${right_file%.ftm10}.fq minlen=25
#    fi
#done
#echo "Done remove read less than 25bp"

#histogram generation
# for index in ${!array[*]}; do
	#bbduk.sh in1=quality_trim_and_filtered_${array[$index]} in2=quality_trim_and_filtered_${array2[$index]} bhist=bhist.txt qhist=qhist.txt gchist=gchist.txt aqhist=aqhist.txt lhist=lhist.txt gcbins=auto
# done
#array6=($(ls $G/*_1.fastq.ftmclean.adpclean.trim.trim.filter))
#for i in ${array6[@]}; do 
#	bbduk.sh in1=${i} in2=$(echo ${i}|sed s/_1/_2/) bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
# 	#All histogram output contents are combined into one file
# 	echo "STOP" $(date)
#	echo ${i} > ${i}.hist.all 
#	echo "bhist" >> ${i}.hist.all
#	cat ${i}.b.hist >> ${i}.hist.all
#	echo "qhist" >> ${i}.hist.all
#   cat ${i}.q.hist	>> ${i}.hist.all
#	echo "gchist" >> ${i}.hist.all
#    cat ${i}.gc.hist >> ${i}.hist.all
#	echo "lhist" >> ${i}.hist.all
#	cat ${i}.l.hist >> ${i}.hist.all
# done

	#lhist = output a read length histogram
	#qhist = per base average quality
	#bhist = output a per-base composition histogram 
	#gchist = output a gc content histogram



#tophat
	#bowtie2-build
#cufflinks

#RSEM

echo "STOP $(date)"
