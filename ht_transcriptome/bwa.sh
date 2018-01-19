#!/bin/bash
#PBS -l walltime=1000:00:00
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -o bwa_out

set -e
echo "START" $(date)

module load BWA/0.7.15-foss-2016b
#bwa code: 
#

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
bac_genomes=${project_home}/bac_genomes
bwa_index=${bac_genomes}/bwa_index/allgenomes_bac_bwa.fna
bwa=${project_home}/bwa
trans2017=${bwa}/trans2017
qc=$project_home/qc/trans2017

cd $trans2017
source $project_home/github/oyster/lib/slack.sh

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

#bwa mem -t 4 -a /net/fs03/data3/marine_diseases_lab/tejashree/Bio_project_SRA/bac_genomes/bwa_index/allgenomes_bac_bwa.fna SRR5357619_1.fastq SRR5357619_2.fastq > SRR5357619.fastq.sam

# ------------------------------------------------------------------------------------------
#bwa align paired end reads to reference genomes of all bacteria 
array1=($(ls -1 ${qc}/${FILES[$PBS_ARRAYID]}_1.fq))
if [ "${#array1[@]}" == "0" ]; then
    echo "($PBS_ARRAYID) ERROR: No input files matching '${qc}/${FILES[$PBS_ARRAYID]}_1.fq'"
    post_slack_message cluster-jobs "ERROR: No input files matching '${qc}/${FILES[$PBS_ARRAYID]}_1.fq'" "$me"
    exit 1
fi

echo "Starting bwa alignment on ${array1[@]}"
post_slack_message cluster-jobs "Starting bwa alignment on ${array1[@]}" "$me"
for left_file in ${array1[@]}; do
    right_file=$(echo ${left_file}|sed s/_1/_2/)
    file=$(basename ${left_file%_1.fq})

    echo "starting bwa mem -t4 -a ${bwa_index} ${left_file} ${right_file} > ${trans2017}/${file}.sam"
    post_slack_message cluster-jobs "bwa mem -t4 -a ${bwa_index} ${left_file} ${right_file} > ${trans2017}/${file}.sam" "${me}.${file}"
    if [ "$debug" ] ; then
        touch ${trans2017}/${file}.sam
    else
        bwa mem -t4 -a ${bwa_index} ${left_file} ${right_file} > ${trans2017}/${file}.sam
    fi
    echo done
    post_slack_message cluster-jobs "done" "${me}.${file}"
done
echo "Done bwa alignment"
post_slack_message cluster-jobs "Done bwa alignment" "$me"
echo "STOP" $(date)
