#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
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

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
cvir_genome=${project_home}/cvir_genome
hisat2=${project_home}/hisat2
qc=${project_home}/qc

cd $hisat2
mkdir -p cvir

#cd /data3/marine_diseases_lab/tejashree/Bio_project_SRA/hisat2/cvir_genome
#cd /data3/marine_diseases_lab/tejashree/Bio_project_SRA/hisat2/trim5
#F=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/PE_fastq/ftm_trim
#S=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/hisat2/cvir_genome

# hisat2-build -f /data3/marine_diseases_lab/tejashree/Bio_project_SRA/hisat2/cvir_genome/cvir.fa  cvir

# -f indicates that the reference input files are FASTA files

#Stay in the directory created in the previous step

#Aligning paired end reads
array1=($(ls -1 ${qc}/*_1.fq))
echo "Starting hisat2 alignment"
for left_file in ${array1[@]}; do
    right_file=$(echo ${left_file}|sed s/_1/_2/)
    file=$(basename ${left_file%_1.fq}) 
    echo starting hisat2 --dta -x ${cvir_genome}/cvir -1 ${left_file} -2 ${right_file} -S ${hisat2}/cvir/${file}.sam
    if [ "$debug" ] ; then
        touch ${hisat2}/cvir/${file}.sam
    else
        hisat2 --dta -x ${cvir_genome}/cvir -1 ${left_file} -2 ${right_file} -S ${hisat2}/cvir/${file}.sam
    fi
    echo done
done
echo "Done hisat2 alignment"
 	#don't need -f because the reads are fastq
	# put -x before the index
	# --dta : Report alignments tailored for transcript assemblers including StringTie.
	 #With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. 
	 #This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

#Aligning single end reads, removed unprocessed PE reads first before proceeding due to same .fq suffix
#array2=($(ls $F/*.fq))

#for i in ${array2[@]}; do
#	hisat2 --dta -x /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/genome_index -U ${i} -S ${i}.sam
#	echo "${i}_DONE"
#done
	#put -x right before the genome index base, and the -U right before unpaired reads

#SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie
array3=($(ls -1 ${hisat2}/cvir/*.sam))
echo "Starting samtools convert"
for i in ${array3[@]}; do
	echo starting samtools sort -o ${i%.sam}.bam ${i}
    if [ "$debug" ] ; then
        touch ${i%.sam}.bam
    else
        samtools sort -o ${i%.sam}.bam ${i}
    fi
    echo done
done
echo "Done samtools convert"
echo "STOP" $(date)
