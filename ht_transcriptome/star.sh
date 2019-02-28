#!/bin/bash
#PBS -l walltime=1000:00:00
#PBS -j oe 
#PBS -o star_out

#To run: 
#qsub -N nameofjob -l nodes=1:ppn=10 -t 1-5 star.sh 

set -e
echo "START" $(date)

module load STAR/2.5.3a-foss-2016b

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA
sra=$project_home/sra
star_genome=$project_home/star/genome/
outdir=$star_genome/test_parent/

FILES=( _ SRR5357620 SRR5357623 SRR5357624 SRR5357625 SRR5357626 )

source $project_home/github/oyster/lib/slack.sh

if [ "$PBS_ARRAYID" ] && [ "$PBS_ARRAYID" -lt "${#FILES[@]}" ] && [ "$PBS_ARRAYID" -ge 0 ] ; then
    me="${PBS_JOBNAME}@${HOSTNAME}-${FILES[$PBS_ARRAYID]}"
    file_prefix="${FILES[$PBS_ARRAYID]}"
    left_file="$sra/${file_prefix}_1.fastq" right_file="$sra/${file_prefix}_2.fastq"

    if [ -f "$left_file" ] && [ -f "$right_file" ] ; then #-f checks if that file exists
        post_slack_message cluster-jobs "$me: Starting STAR job 'STAR outFileNamePrefix:$outdir/$file_prefix outdir:$outdir"

        STAR --runThreadN 12 --genomeDir $outdir --sjdbGTFtagExonParentTranscript Parent \
             --outFileNamePrefix $outdir/$file_prefix \
             --sjdbGTFfile $star_genome/ref_C_virginica-3.0_top_level.gff3 --sjdbOverhang 100 \
             --readFilesIn $left_file $right_file

        post_slack_message cluster-jobs "$me: Done STAR job 'STAR outFileNamePrefix:$outdir/$file_prefix outdir:$outdir"
    else
        post_slack_message cluster-jobs "$me: '$left_file' '$right_file' not found"    
    fi
fi

#STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data3/marine_diseases_lab/tejashree/Bio_project_SRA/star/genome/test_parent --sjdbGTFtagExonParentTranscript Parent --genomeFastaFiles /data3/marine_diseases_lab/tejashree/Bio_project_SRA/star/genome/cvir_edited.fa 



#STAR --runThreadN 12 --genomeDir /data3/marine_diseases_lab/tejashree/Bio_project_SRA/star/genome/ --sjdbGTFfile /data3/marine_diseases_lab/tejashree/Bio_project_SRA/cvir_genome/gffread.gtf --sjdbOverhang 100 --readFilesIn /data3/marine_diseases_lab/tejashree/Bio_project_SRA/sra/SRR5357619_1.fastq /data3/marine_diseases_lab/tejashree/Bio_project_SRA/sra/SRR5357619_2.fastq


echo "done"
